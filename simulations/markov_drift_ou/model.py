"""Core model: OU-drift of transition-matrix logits (OU class, § 5.2).

Implements:
- Continuous OU-drift on K = k(k-1) off-diagonal logits theta_{ij}(t) of a
  k-state row-stochastic Markov transition matrix.
- Euler-Maruyama integration with Delta t = 1.
- Softmax row-restoration (diagonal logit fixed at 0, off-diagonal free).
- Online maximum-likelihood Bayesian learner with decaying step eta_t = c/sqrt(t)
  (approximates the Bayes posterior under OU prior; see § 5.2 derivation).
- Predictive information functionals identical to markov_drift/.
- Energy budget with stationary R = 1 plus storage cost nu(t)*|M|/tau_{1/2}.

All entropies in nats. Energetic cost normalized so k_B T ln 2 = 1.
"""

from __future__ import annotations

import numpy as np


# ---------- OU on logits ----------

def init_logits(k: int, sigma0: float, rng: np.random.Generator) -> np.ndarray:
    """Initial logits theta_{ij}(0) ~ N(0, sigma0^2) for i != j; diagonal = 0.

    Returns (k, k) array; entry (i, i) is held at 0.0 by convention.
    """
    theta = rng.normal(0.0, sigma0, size=(k, k))
    np.fill_diagonal(theta, 0.0)
    return theta


def softmax_rows_full(theta: np.ndarray) -> np.ndarray:
    """Row-wise softmax over all k entries (including diagonal = self-loop).

    This produces a row-stochastic P with self-transition prob P_{ii} controlled
    by theta_{ii} (held at 0).
    """
    # numerical-safe softmax
    m = theta.max(axis=1, keepdims=True)
    e = np.exp(theta - m)
    return e / e.sum(axis=1, keepdims=True)


def ou_step(theta: np.ndarray, theta_star: np.ndarray, lam: float, sigma: float,
            rng: np.random.Generator, dt: float = 1.0) -> np.ndarray:
    """Euler-Maruyama update of off-diagonal OU logits.

    theta_{ij}(t+dt) = theta_{ij}(t) - lam*dt*(theta - theta*) + sigma*sqrt(dt)*xi
    Diagonal entries held at theta_star value (= 0 by default).
    """
    k = theta.shape[0]
    noise = rng.normal(0.0, 1.0, size=(k, k))
    drift = -lam * dt * (theta - theta_star)
    new = theta + drift + sigma * np.sqrt(dt) * noise
    # keep diagonal pinned to theta_star (the diagonal logit anchor)
    di = np.diag_indices(k)
    new[di] = theta_star[di]
    return new


# ---------- Predictive information ----------

def stationary_distribution(P: np.ndarray) -> np.ndarray:
    k = P.shape[0]
    try:
        vals, vecs = np.linalg.eig(P.T)
        idx = int(np.argmin(np.abs(vals - 1.0)))
        v = np.real(vecs[:, idx])
        v = np.abs(v)
        s = v.sum()
        if s > 0 and np.isfinite(s):
            return v / s
    except np.linalg.LinAlgError:
        pass
    return np.full(k, 1.0 / k)


def entropy_of_distribution(p: np.ndarray, eps: float = 1e-12) -> float:
    p = np.clip(p, eps, 1.0)
    return float(-np.sum(p * np.log(p)))


def kl_div_rows(P: np.ndarray, Q: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    P = np.clip(P, eps, 1.0)
    Q = np.clip(Q, eps, 1.0)
    return np.sum(P * (np.log(P) - np.log(Q)), axis=1)


def matrix_power(P: np.ndarray, tau: int) -> np.ndarray:
    if tau < 1:
        raise ValueError("tau must be >= 1")
    result = np.eye(P.shape[0])
    base = P.copy()
    n = tau
    while n > 0:
        if n & 1:
            result = result @ base
        n >>= 1
        if n:
            base = base @ base
    return result


def optimal_predictive_info(P: np.ndarray, tau: int) -> float:
    pi = stationary_distribution(P)
    Ptau = matrix_power(P, tau)
    H_marginal = entropy_of_distribution(pi)
    H_cond = float(np.sum(pi * np.array([entropy_of_distribution(Ptau[i]) for i in range(P.shape[0])])))
    return max(H_marginal - H_cond, 0.0)


def observer_predictive_info(P_true: np.ndarray, P_hat: np.ndarray, tau: int) -> float:
    """I_obs = I_opt(P_true, tau) - E_{x ~ pi_true}[KL(P_true^tau || P_hat^tau)]."""
    pi_true = stationary_distribution(P_true)
    Pt_tau = matrix_power(P_true, tau)
    Ph_tau = matrix_power(P_hat, tau)
    I_opt = optimal_predictive_info(P_true, tau)
    kl = float(np.sum(pi_true * kl_div_rows(Pt_tau, Ph_tau)))
    return max(I_opt - kl, 0.0)


# ---------- Online Bayesian learner ----------

class OnlineMLLogitEstimator:
    """Online maximum-likelihood estimator of transition-matrix logits.

    Approximates the Bayes posterior under an OU prior with stationary variance
    sigma^2/(2 lambda) by gradient ascent on the conditional log-likelihood with
    decaying step eta_t = c / sqrt(t).

    Gradient of log P_{i,.}(x_{t+1} | theta) w.r.t. theta_{ij} when x_t = i:
        grad_{ij} = 1[x_{t+1} = j] - P_{ij}(theta).
    For rows other than i the gradient is 0 (observation provides no info on them
    at this step). Diagonal entry theta_{ii} held fixed = 0.

    Returns hat_theta and hat_P via softmax_rows_full(hat_theta).
    """

    def __init__(self, k: int, c: float = 1.0, theta0: float = 0.0):
        self.k = k
        self.c = c
        self.theta = np.full((k, k), theta0, dtype=np.float64)
        np.fill_diagonal(self.theta, 0.0)
        self._step = 0

    def update(self, x_prev: int, x_next: int) -> None:
        self._step += 1
        eta = self.c / np.sqrt(self._step)
        # Current row probabilities under softmax over all k entries
        row = self.theta[x_prev]
        m = row.max()
        e = np.exp(row - m)
        p = e / e.sum()
        grad = -p
        grad[x_next] += 1.0
        # do not update the diagonal entry (kept as anchor)
        grad[x_prev] = 0.0
        self.theta[x_prev] = row + eta * grad

    def estimate_P(self) -> np.ndarray:
        return softmax_rows_full(self.theta)

    def estimate_theta(self) -> np.ndarray:
        return self.theta.copy()


# ---------- Simulation driver ----------

def simulate_run(
    k: int,
    T: int,
    lam: float,
    sigma: float,
    tau_pred: int,
    rng: np.random.Generator,
    measure_every: int = 100,
    theta0_std: float = 0.1,
    R: float = 1.0,
    N0: float = 0.0,
    tau_half: float = 100.0,
    learn_c: float = 1.0,
) -> dict:
    """Run a single OU-drift simulation.

    Returns a dict with arrays sampled every `measure_every` steps:
      times, I_pred, I_opt, nu, N_max_R, N_max_full, eta_L_R, eta_L_full,
      theta_norm (||theta(t) - hat_theta(t)||_F for diagnostic).

    Two N_max series:
      - N_max_R(t)    = R * t                          (matches theory derivation)
      - N_max_full(t) = sum_s [R + nu(s) * K / tau_half] (full spec budget)
    Two efficiencies analogously.
    """
    K_params = k * (k - 1)

    # initialise true OU logits and ideal-learner estimator
    theta_star = np.zeros((k, k))
    theta = init_logits(k, theta0_std, rng)
    estimator = OnlineMLLogitEstimator(k, c=learn_c)

    # initial state: stationary of current P
    P_true = softmax_rows_full(theta)
    pi0 = stationary_distribution(P_true)
    x = int(rng.choice(k, p=pi0))

    n_samples = T // measure_every
    times = np.zeros(n_samples)
    I_pred_arr = np.zeros(n_samples)
    I_opt_arr = np.zeros(n_samples)
    nu_arr = np.zeros(n_samples)
    NmaxR_arr = np.zeros(n_samples)
    Nmaxfull_arr = np.zeros(n_samples)
    etaR_arr = np.zeros(n_samples)
    etafull_arr = np.zeros(n_samples)
    theta_err_arr = np.zeros(n_samples)

    # incremental budget accumulators
    Nmax_full = float(N0)
    rate_store_coef = K_params / tau_half  # multiplied by nu(t) each step
    last_nu_for_step = 0.0

    sample_idx = 0
    for t in range(1, T + 1):
        # advance OU logits
        theta = ou_step(theta, theta_star, lam, sigma, rng, dt=1.0)
        P_true = softmax_rows_full(theta)

        # advance Markov chain
        x_next = int(rng.choice(k, p=P_true[x]))
        estimator.update(x, x_next)
        x = x_next

        # accumulate full N_max using last measured nu (piecewise-constant
        # between measurements; for t < first measurement use 0)
        Nmax_full += R + last_nu_for_step * rate_store_coef

        # measurement
        if t % measure_every == 0 and sample_idx < n_samples:
            P_hat = estimator.estimate_P()
            I_opt = optimal_predictive_info(P_true, tau_pred)
            I_obs = observer_predictive_info(P_true, P_hat, tau_pred)
            nu = 1.0 - (I_obs / I_opt) if I_opt > 1e-9 else 0.0
            NmaxR = N0 + R * t
            theta_hat = estimator.estimate_theta()
            theta_err = float(np.linalg.norm(theta - theta_hat))

            times[sample_idx] = t
            I_pred_arr[sample_idx] = I_obs
            I_opt_arr[sample_idx] = I_opt
            nu_arr[sample_idx] = nu
            NmaxR_arr[sample_idx] = NmaxR
            Nmaxfull_arr[sample_idx] = Nmax_full
            etaR_arr[sample_idx] = I_obs / NmaxR if NmaxR > 0 else 0.0
            etafull_arr[sample_idx] = I_obs / Nmax_full if Nmax_full > 0 else 0.0
            theta_err_arr[sample_idx] = theta_err
            last_nu_for_step = float(nu)
            sample_idx += 1

    return dict(
        times=times[:sample_idx],
        I_pred=I_pred_arr[:sample_idx],
        I_opt=I_opt_arr[:sample_idx],
        nu=nu_arr[:sample_idx],
        N_max_R=NmaxR_arr[:sample_idx],
        N_max_full=Nmaxfull_arr[:sample_idx],
        eta_L_R=etaR_arr[:sample_idx],
        eta_L_full=etafull_arr[:sample_idx],
        theta_err=theta_err_arr[:sample_idx],
    )


def average_runs(runs: list[dict]) -> dict:
    keys = ["I_pred", "I_opt", "nu", "N_max_R", "N_max_full",
            "eta_L_R", "eta_L_full", "theta_err"]
    times = runs[0]["times"]
    out = {"times": times}
    for key in keys:
        stacked = np.stack([r[key] for r in runs], axis=0)
        out[key] = stacked.mean(axis=0)
        out[key + "_std"] = stacked.std(axis=0)
    return out
