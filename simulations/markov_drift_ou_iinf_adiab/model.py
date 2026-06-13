"""Vectorised model for iter-006 Step 2: OU + I_infinity adiabatic scan.

Same physics as ``../markov_drift_ou_iinf/model.py`` (OU-drift k=8 Markov chain
with K = k(k-1) = 56 off-diagonal logits, ideal online ML learner, cumulative
L_excess(t) BNT redundancy measurement) but with all N_RUNS realisations packed
into leading-axis arrays of shape (N, k, k) so the per-step inner loop runs in
vectorised numpy. This is necessary for T = 2e5 with N ~ 20 realisations to
finish in minutes rather than hours.

The estimator state is a (N, k, k) logit tensor; the chain state is an int
array of length N; the OU drift is an Euler-Maruyama update broadcast over all
realisations; sampling from each chain's current row of P uses inverse-CDF on
the per-row cumulative distribution.

Outputs (per measurement time) are arrays of shape (N, n_samples).
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# Vectorised primitives                                                       #
# --------------------------------------------------------------------------- #

def init_logits_batched(N: int, k: int, sigma0: float,
                        rng: np.random.Generator) -> np.ndarray:
    """Initial logits, shape (N, k, k). Diagonal pinned at 0."""
    theta = rng.normal(0.0, sigma0, size=(N, k, k))
    # zero diagonals across the batch
    idx = np.arange(k)
    theta[:, idx, idx] = 0.0
    return theta


def softmax_rows_batched(theta: np.ndarray) -> np.ndarray:
    """Row-wise softmax, shape (N, k, k) -> (N, k, k)."""
    m = theta.max(axis=-1, keepdims=True)
    e = np.exp(theta - m)
    return e / e.sum(axis=-1, keepdims=True)


def ou_step_batched(theta: np.ndarray, lam: float, sigma: float,
                    rng: np.random.Generator, dt: float = 1.0) -> np.ndarray:
    """Euler-Maruyama on the batched OU logits.

    Theta* = 0 (matches Step 1). Diagonal pinned to 0 after the update.
    """
    N, k, _ = theta.shape
    noise = rng.standard_normal(size=(N, k, k))
    new = theta - lam * dt * theta + sigma * np.sqrt(dt) * noise
    idx = np.arange(k)
    new[:, idx, idx] = 0.0
    return new


def sample_next_states(P: np.ndarray, x: np.ndarray,
                       rng: np.random.Generator) -> np.ndarray:
    """Sample x_next ~ P[n, x[n], :] for each realisation n.

    P: (N, k, k); x: (N,) int. Returns (N,) int.
    """
    N, k, _ = P.shape
    # Pick out the current row for each realisation, shape (N, k)
    rows = P[np.arange(N), x, :]
    # Inverse-CDF sampling: u ~ U[0,1), find first j s.t. cumsum >= u
    cdf = np.cumsum(rows, axis=1)
    # numerical safety: snap last column to 1.0 + eps
    cdf[:, -1] = np.maximum(cdf[:, -1], 1.0)
    u = rng.random(size=(N, 1))
    x_next = (cdf >= u).argmax(axis=1)
    return x_next.astype(np.int64)


def online_ml_update_batched(theta: np.ndarray, x_prev: np.ndarray,
                             x_next: np.ndarray, step_count: int,
                             c: float) -> None:
    """In-place online ML update of batched logits.

    Per Step 1 estimator:
      grad_{ij} = 1[x_next == j] - P_{ij}    for row i = x_prev
      diagonal grad pinned to 0
      step size eta = c / sqrt(step_count)
    """
    N, k, _ = theta.shape
    eta = c / np.sqrt(step_count)
    rows = theta[np.arange(N), x_prev, :]  # (N, k)
    m = rows.max(axis=1, keepdims=True)
    e = np.exp(rows - m)
    p = e / e.sum(axis=1, keepdims=True)   # (N, k)
    grad = -p
    grad[np.arange(N), x_next] += 1.0
    # zero the diagonal entry per realisation (the (x_prev[n], x_prev[n]) anchor)
    grad[np.arange(N), x_prev] = 0.0
    theta[np.arange(N), x_prev, :] = rows + eta * grad


def per_row_log_density(P: np.ndarray, x_prev: np.ndarray, x_next: np.ndarray,
                        eps: float = 1e-12) -> np.ndarray:
    """log P[n, x_prev[n], x_next[n]] for each realisation n. Returns (N,)."""
    N = P.shape[0]
    vals = P[np.arange(N), x_prev, x_next]
    return np.log(np.maximum(vals, eps))


# --------------------------------------------------------------------------- #
# Stationary distribution (per realisation, batched eigendecomp at sample time)
# --------------------------------------------------------------------------- #

def stationary_distribution_batched(P: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Stationary pi for each P[n], using a few iterations of power method
    starting from uniform. The chain mixes fast (lambda_2 small for typical
    softmax logits) so 200 iterations is more than enough for our k=8.
    """
    N, k, _ = P.shape
    pi = np.full((N, k), 1.0 / k)
    for _ in range(200):
        # batched (N,k) <- (N,k) * (N,k,k):  use einsum to enforce broadcast over N
        pi = np.einsum("ni,nij->nj", pi, P)
        pi = pi / pi.sum(axis=1, keepdims=True)
    return pi


def optimal_predictive_info_batched(P: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """I_opt = H(pi) - sum_i pi_i H(P[i,:]), one value per realisation, shape (N,)."""
    pi = stationary_distribution_batched(P)  # (N, k)
    H_marg = -np.sum(pi * np.log(np.maximum(pi, eps)), axis=1)  # (N,)
    H_row = -np.sum(P * np.log(np.maximum(P, eps)), axis=2)     # (N, k)
    H_cond = np.sum(pi * H_row, axis=1)                          # (N,)
    return np.maximum(H_marg - H_cond, 0.0)


def observer_predictive_info_batched(P_true: np.ndarray, P_hat: np.ndarray,
                                     eps: float = 1e-12) -> np.ndarray:
    """I_obs = I_opt(P_true) - E_pi[KL(P_true[i,:] || P_hat[i,:])], shape (N,)."""
    pi = stationary_distribution_batched(P_true)            # (N, k)
    P_true_c = np.maximum(P_true, eps)
    P_hat_c = np.maximum(P_hat, eps)
    kl_rows = np.sum(P_true_c * (np.log(P_true_c) - np.log(P_hat_c)), axis=2)  # (N, k)
    kl = np.sum(pi * kl_rows, axis=1)                         # (N,)
    I_opt = optimal_predictive_info_batched(P_true)
    return np.maximum(I_opt - kl, 0.0)


# --------------------------------------------------------------------------- #
# Top-level simulator                                                         #
# --------------------------------------------------------------------------- #

def simulate_batch_iinf(
    k: int,
    T: int,
    lam: float,
    sigma: float,
    n_runs: int,
    rng: np.random.Generator,
    measure_every: int = 100,
    theta0_std: float = 0.1,
    learn_c: float = 3.0,
    R: float = 1.0,
) -> dict:
    """Vectorised version of simulate_run_iinf for n_runs realisations.

    Returns dict with arrays:
      times                shape (n_samples,)
      L_excess_cum_mean    shape (n_samples,)    mean over realisations
      L_excess_cum_std     shape (n_samples,)
      G_obs_cum_mean       shape (n_samples,)
      cumulative_I_opt_mean shape (n_samples,)
      I_obs_1step_mean     shape (n_samples,)
      I_opt_1step_mean     shape (n_samples,)
      nu_mean              shape (n_samples,)
      theta_err_mean       shape (n_samples,)
      n_runs               int
    """
    K_params = k * (k - 1)
    N = n_runs

    # initial logits and learner state
    theta = init_logits_batched(N, k, theta0_std, rng)
    theta_hat = np.zeros((N, k, k))  # learner starts from uniform-row logits

    # initial state x[n] ~ pi for each realisation
    P_true = softmax_rows_batched(theta)
    pi0 = stationary_distribution_batched(P_true)
    # Sample by inverse-CDF from pi0 per realisation
    cdf0 = np.cumsum(pi0, axis=1)
    cdf0[:, -1] = np.maximum(cdf0[:, -1], 1.0)
    u0 = rng.random(size=(N, 1))
    x = (cdf0 >= u0).argmax(axis=1).astype(np.int64)

    n_samples = T // measure_every
    times = np.zeros(n_samples)
    L_excess_cum_arr = np.zeros((N, n_samples))
    G_obs_cum_arr = np.zeros((N, n_samples))
    cum_I_opt_arr = np.zeros((N, n_samples))
    I_obs_arr = np.zeros((N, n_samples))
    I_opt_arr = np.zeros((N, n_samples))
    nu_arr = np.zeros((N, n_samples))
    theta_err_arr = np.zeros((N, n_samples))

    L_excess_cum = np.zeros(N)
    G_obs_cum = np.zeros(N)
    cum_I_opt = np.zeros(N)
    log_k = float(np.log(k))

    sample_idx = 0
    for t in range(1, T + 1):
        # advance OU logits
        theta = ou_step_batched(theta, lam, sigma, rng, dt=1.0)
        P_true = softmax_rows_batched(theta)

        # advance Markov chain (per-realisation)
        x_next = sample_next_states(P_true, x, rng)

        # PRE-UPDATE learner P_hat: same convention as Step 1
        P_hat = softmax_rows_batched(theta_hat)
        log_p_true = per_row_log_density(P_true, x, x_next)
        log_p_hat = per_row_log_density(P_hat, x, x_next)
        L_excess_cum += (log_p_true - log_p_hat)
        G_obs_cum += (log_p_hat + log_k)

        # accumulate instantaneous I_opt (only sample at measurement times to save cost)
        # learner online update
        online_ml_update_batched(theta_hat, x, x_next, t, learn_c)
        x = x_next

        if t % measure_every == 0 and sample_idx < n_samples:
            # measurement: compute I_opt, I_obs, nu, theta_err
            I_opt = optimal_predictive_info_batched(P_true)
            # we accumulate I_opt only at sample times, scaled by measure_every
            # so cum_I_opt approximates sum_s I_opt_s at coarser resolution
            cum_I_opt += measure_every * I_opt
            P_hat_post = softmax_rows_batched(theta_hat)
            I_obs = observer_predictive_info_batched(P_true, P_hat_post)
            with np.errstate(divide="ignore", invalid="ignore"):
                nu = np.where(I_opt > 1e-9, 1.0 - I_obs / I_opt, 0.0)
            theta_err = np.linalg.norm(theta - theta_hat, axis=(1, 2))

            times[sample_idx] = t
            L_excess_cum_arr[:, sample_idx] = L_excess_cum
            G_obs_cum_arr[:, sample_idx] = G_obs_cum
            cum_I_opt_arr[:, sample_idx] = cum_I_opt
            I_obs_arr[:, sample_idx] = I_obs
            I_opt_arr[:, sample_idx] = I_opt
            nu_arr[:, sample_idx] = nu
            theta_err_arr[:, sample_idx] = theta_err
            sample_idx += 1

    times = times[:sample_idx]
    L_excess_cum_arr = L_excess_cum_arr[:, :sample_idx]
    G_obs_cum_arr = G_obs_cum_arr[:, :sample_idx]
    cum_I_opt_arr = cum_I_opt_arr[:, :sample_idx]
    I_obs_arr = I_obs_arr[:, :sample_idx]
    I_opt_arr = I_opt_arr[:, :sample_idx]
    nu_arr = nu_arr[:, :sample_idx]
    theta_err_arr = theta_err_arr[:, :sample_idx]

    return dict(
        n_runs=N,
        k=k,
        K_params=K_params,
        times=times,
        L_excess_cum=L_excess_cum_arr.mean(axis=0),
        L_excess_cum_std=L_excess_cum_arr.std(axis=0),
        L_excess_cum_full=L_excess_cum_arr,   # keep raw for per-run fits
        G_obs_cum=G_obs_cum_arr.mean(axis=0),
        cumulative_I_opt=cum_I_opt_arr.mean(axis=0),
        I_obs_1step=I_obs_arr.mean(axis=0),
        I_opt_1step=I_opt_arr.mean(axis=0),
        nu=nu_arr.mean(axis=0),
        nu_std=nu_arr.std(axis=0),
        theta_err=theta_err_arr.mean(axis=0),
    )
