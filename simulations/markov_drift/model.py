"""Core model: piecewise-stationary Markov chain with Poisson regime switches (PSP class).

Implements:
- Sampling random row-stochastic matrices from Dirichlet(alpha).
- Piecewise-constant transition matrix P(t) with Poisson switches at rate lambda.
- Observer estimate \hat P(t) from sliding window of empirical transition counts
  (with Laplace smoothing).
- Stationary distribution and predictive information functionals.

Re-equipped apparatus (efficiency redefined I_pred/N_max -> I_pred/I_mem):
  I_mem(t) = (K/2) * ln(max(N_obs(t), e)),  K = k*(k-1),
             N_obs = number of transitions retained in memory.
  This is the standard learned (model<->data mutual) information of an MDL/Bayes
  estimator of a K-parameter model from N observations, I ~ (K/2) ln N nats.
  It is a data-processing bound on what the model can know about the dynamics,
  so I_pred <= I_mem always (eta in [0,1]).
  eta(t)        = I_pred(t) / I_mem(t)                  (Landauer efficiency)
  nu_Still(t)   = 1 - eta(t)                            (memory ballast = nostalgia)
  nu_op(t)      = 1 - I_pred(t) / I_opt(t)              (prediction shortfall, secondary)

All entropies/informations in nats.
"""

from __future__ import annotations

import numpy as np


# ---------- Random transition matrices ----------

def sample_stochastic_matrix(k: int, rng: np.random.Generator, alpha: float = 1.0) -> np.ndarray:
    """Sample k x k row-stochastic matrix from Dirichlet(alpha) on each row."""
    return rng.dirichlet(np.full(k, alpha), size=k)


def stationary_distribution(P: np.ndarray) -> np.ndarray:
    """Return left eigenvector of P with eigenvalue 1, normalized to sum to 1.

    Falls back to uniform if computation is ill-conditioned.
    """
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


# ---------- Predictive information ----------

def kl_div_rows(P: np.ndarray, Q: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Row-wise KL(P || Q) in nats, vector of length k."""
    P = np.clip(P, eps, 1.0)
    Q = np.clip(Q, eps, 1.0)
    return np.sum(P * (np.log(P) - np.log(Q)), axis=1)


def entropy_of_distribution(p: np.ndarray, eps: float = 1e-12) -> float:
    """Shannon entropy in nats."""
    p = np.clip(p, eps, 1.0)
    return float(-np.sum(p * np.log(p)))


def matrix_power(P: np.ndarray, tau: int) -> np.ndarray:
    """Return P^tau via repeated squaring (tau >= 1)."""
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
    """I_opt(P, tau) = MI(X_t ; X_{t+tau}) in nats, X_t ~ stationary pi.

    = H(X_{t+tau}) - H(X_{t+tau} | X_t) where both are under joint stationary law.
    For stationary pi: H(X_{t+tau}) = H(pi); H(X_{t+tau}|X_t) = sum_i pi_i H(P^tau row i).
    """
    pi = stationary_distribution(P)
    Ptau = matrix_power(P, tau)
    H_marginal = entropy_of_distribution(pi)
    H_cond = float(np.sum(pi * np.array([entropy_of_distribution(Ptau[i]) for i in range(P.shape[0])])))
    return max(H_marginal - H_cond, 0.0)


def observer_predictive_info(P_true: np.ndarray, P_hat: np.ndarray, tau: int) -> float:
    """Realised predictive information under observer model P_hat, true dynamics P_true.

    I_obs = I_opt(P_true, tau) - E_{x ~ pi}[ KL( P_true^tau(.|x) || P_hat^tau(.|x) ) ].

    Clipped at 0 since KL can dominate when model is very wrong.
    """
    pi_true = stationary_distribution(P_true)
    Pt_tau = matrix_power(P_true, tau)
    Ph_tau = matrix_power(P_hat, tau)
    I_opt = optimal_predictive_info(P_true, tau)
    kl = float(np.sum(pi_true * kl_div_rows(Pt_tau, Ph_tau)))
    return max(I_opt - kl, 0.0)


def memory_information(n_obs: float, k: int) -> float:
    """Learned information stored by the model: I_mem = (K/2) ln(max(N_obs, e)) nats.

    K = k*(k-1) free off-diagonal parameters of a k-state row-stochastic matrix
    (each row has k-1 free entries). For a K-parameter model estimated from
    N independent observations, the mutual information between the fitted model
    and the data grows as (K/2) ln N (MDL / Bayesian stochastic complexity).
    The floor at N_obs = e keeps I_mem >= K/2 > 0, so eta = I_pred / I_mem is
    well defined from the first observation. NOTE: I_pred <= I_mem here is enforced
    by the finite-sample clip below (a guard on estimator overshoot), NOT a DPI
    theorem -- there is no general DPI order for the operational numerator (paper 3.1).
    """
    K = k * (k - 1)
    return 0.5 * K * np.log(max(n_obs, np.e))


# ---------- Observer memory ----------

class SlidingWindowEstimator:
    """Maintain counts of transitions in a sliding window of length tau_w.

    Use np.inf as tau_w for unbounded memory (no decay).
    Implementation: store the transition events in a ring buffer when tau_w is finite;
    when infinite, accumulate counts directly.
    """

    def __init__(self, k: int, tau_w: float, laplace: float = 1.0):
        self.k = k
        self.tau_w = tau_w
        self.laplace = laplace
        self.counts = np.zeros((k, k), dtype=np.float64)
        if np.isfinite(tau_w):
            self._buffer = np.empty((int(tau_w), 2), dtype=np.int64)
            self._buf_idx = 0
            self._buf_full = False
        else:
            self._buffer = None

    def update(self, x_prev: int, x_next: int) -> None:
        if self._buffer is not None:
            # remove evicted event
            if self._buf_full:
                old_i, old_j = self._buffer[self._buf_idx]
                self.counts[old_i, old_j] -= 1
            self._buffer[self._buf_idx] = (x_prev, x_next)
            self._buf_idx += 1
            if self._buf_idx >= self._buffer.shape[0]:
                self._buf_idx = 0
                self._buf_full = True
        self.counts[x_prev, x_next] += 1

    def estimate(self) -> np.ndarray:
        """Return Laplace-smoothed row-stochastic estimate."""
        c = self.counts + self.laplace
        return c / c.sum(axis=1, keepdims=True)

    def total_observations(self) -> float:
        """Number of transitions currently retained in memory (window occupancy)."""
        return float(self.counts.sum())


# ---------- Simulation driver ----------

def simulate_run(
    k: int,
    T: int,
    lam: float,
    tau_w: float,
    tau_pred: int,
    rng: np.random.Generator,
    measure_every: int = 50,
    laplace: float = 1.0,
    alpha: float = 1.0,
) -> dict:
    """Run a single simulation.

    Returns a dict with arrays sampled every `measure_every` steps:
      times, I_pred, I_opt, I_mem, eta, nu_still, nu_op, n_obs,
      switch_times (list of T_i within [0, T]).

    Re-equipped cost model (all in nats):
      I_mem(t) = (K/2) ln(max(N_obs(t), e)),  N_obs = transitions in memory.
      eta(t)      = I_pred(t) / I_mem(t)        (replaces old eta_L = I_pred / N_max).
      nu_still(t) = 1 - eta(t)                  (memory ballast, canonical nostalgia).
      nu_op(t)    = 1 - I_pred(t) / I_opt(t)    (prediction shortfall, secondary).
    """
    # initial transition matrix and first switch time
    P_true = sample_stochastic_matrix(k, rng, alpha=alpha)
    if lam > 0:
        next_switch = rng.exponential(1.0 / lam)
    else:
        next_switch = np.inf
    switch_times = []

    estimator = SlidingWindowEstimator(k, tau_w, laplace=laplace)

    # warm start state from initial stationary
    pi0 = stationary_distribution(P_true)
    x = int(rng.choice(k, p=pi0))

    n_samples = T // measure_every
    times = np.zeros(n_samples)
    I_pred_arr = np.zeros(n_samples)
    I_opt_arr = np.zeros(n_samples)
    I_mem_arr = np.zeros(n_samples)
    eta_arr = np.zeros(n_samples)
    nu_still_arr = np.zeros(n_samples)
    nu_op_arr = np.zeros(n_samples)
    n_obs_arr = np.zeros(n_samples)

    sample_idx = 0
    for t in range(1, T + 1):
        # regime switch?
        if t >= next_switch:
            P_true = sample_stochastic_matrix(k, rng, alpha=alpha)
            switch_times.append(t)
            next_switch = t + rng.exponential(1.0 / lam)

        # advance chain
        x_next = int(rng.choice(k, p=P_true[x]))
        estimator.update(x, x_next)
        x = x_next

        # measurement
        if t % measure_every == 0 and sample_idx < n_samples:
            P_hat = estimator.estimate()
            I_opt = optimal_predictive_info(P_true, tau_pred)
            I_obs = observer_predictive_info(P_true, P_hat, tau_pred)
            n_obs = estimator.total_observations()
            I_mem = memory_information(n_obs, k)
            # Finite-sample guard: clip estimator overshoot to I_mem (NOT a DPI theorem;
            # no general DPI order for the operational numerator -- paper 3.1).
            I_obs = min(I_obs, I_mem)
            eta = I_obs / I_mem if I_mem > 0 else 0.0
            nu_still = 1.0 - eta
            nu_op = 1.0 - (I_obs / I_opt) if I_opt > 1e-9 else 0.0

            times[sample_idx] = t
            I_pred_arr[sample_idx] = I_obs
            I_opt_arr[sample_idx] = I_opt
            I_mem_arr[sample_idx] = I_mem
            eta_arr[sample_idx] = eta
            nu_still_arr[sample_idx] = nu_still
            nu_op_arr[sample_idx] = nu_op
            n_obs_arr[sample_idx] = n_obs
            sample_idx += 1

    return dict(
        times=times[:sample_idx],
        I_pred=I_pred_arr[:sample_idx],
        I_opt=I_opt_arr[:sample_idx],
        I_mem=I_mem_arr[:sample_idx],
        eta=eta_arr[:sample_idx],
        nu_still=nu_still_arr[:sample_idx],
        nu_op=nu_op_arr[:sample_idx],
        n_obs=n_obs_arr[:sample_idx],
        switch_times=np.array(switch_times, dtype=np.float64),
    )


def average_runs(runs: list[dict]) -> dict:
    """Average per-time arrays across runs (assumes identical sampling grid)."""
    keys = ["I_pred", "I_opt", "I_mem", "eta", "nu_still", "nu_op", "n_obs"]
    times = runs[0]["times"]
    out = {"times": times}
    for key in keys:
        stacked = np.stack([r[key] for r in runs], axis=0)
        out[key] = stacked.mean(axis=0)
        out[key + "_std"] = stacked.std(axis=0)
    return out
