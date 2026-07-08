"""Core model: OU-drift of transition-matrix logits (OU class, § 6.2).

Implements:
- Continuous OU-drift on K = k(k-1) off-diagonal logits theta_{ij}(t) of a
  k-state row-stochastic Markov transition matrix.
- Euler-Maruyama integration with Delta t = 1.
- Softmax row-restoration (diagonal logit fixed at 0, off-diagonal free).
- TWO online maximum-likelihood learners on the SAME trajectory (the central
  contrast of § 6.2):
    * Robbins-Monro decaying step eta_t = c/sqrt(t) -- FREEZES under drift;
    * constant step eta_t = s -- TRACKS the moving target to a finite floor.
  NO sliding window: each online learner has seen N_obs = t observations at time t.
- Predictive information functionals identical to markov_drift/.

Re-equipped apparatus (efficiency redefined I_pred/N_max -> I_pred/I_mem):
  I_mem(t) = (K/2) * ln(max(N_obs, e)),  K = k*(k-1),  N_obs = t,
             the number of transition observations seen by the online learner.
  This is the standard learned (model<->data mutual) information of an MDL/Bayes
  estimator of a K-parameter model from N observations, I ~ (K/2) ln N nats.
  It is a data-processing bound on what the model can know about the dynamics,
  so I_pred <= I_mem always (eta in [0,1]).
  eta(t)        = I_pred(t) / I_mem(t)                  (Landauer efficiency)
  nu_Still(t)   = 1 - eta(t)                            (memory ballast, secondary)
  nu_op(t)      = 1 - I_pred(t) / I_opt(t)              (operational nostalgia, CANONICAL)

The energy-budget denominator N_max = R*t (and tau_half storage cost, R, N0,
rate_store_coef) has been REMOVED: it was an artefact of the old efficiency
definition. There is no slope K/(2R) to fit any more.

All entropies/informations in nats.
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


def memory_information(n_obs: float, k: int) -> float:
    """Learned information stored by the model: I_mem = (K/2) ln(max(N_obs, e)) nats.

    K = k*(k-1) free off-diagonal parameters of a k-state row-stochastic matrix
    (each row has k-1 free entries). For a K-parameter model estimated from
    N independent observations, the mutual information between the fitted model
    and the data grows as (K/2) ln N (MDL / Bayesian stochastic complexity).
    The floor at N_obs = e keeps I_mem >= K/2 > 0, so eta = I_pred / I_mem is
    well defined from the first observation and I_pred <= I_mem (DPI) holds.

    For the online learner with no sliding window, N_obs(t) = t.
    """
    K = k * (k - 1)
    return 0.5 * K * np.log(max(n_obs, np.e))


# ---------- Online Bayesian learner ----------

class OnlineMLLogitEstimator:
    """Online maximum-likelihood estimator of transition-matrix logits.

    Two step-size schedules are supported on the SAME gradient update; which one
    is used is a pure design choice of the learner (the central contrast of
    § 6.2):

    - Robbins-Monro / decaying step (default, const_step=None):
        eta_t = c / sqrt(t).
      Approximates the Bayes posterior under an OU prior with stationary
      variance sigma^2/(2 lambda). Under a NON-STATIONARY target it FREEZES:
      as t grows the step vanishes, the estimate stops tracking the drift, and
      the operational nostalgia nu_op climbs toward 1.

    - Constant step (const_step=s):
        eta_t = s  (fixed).
      Keeps tracking the moving target indefinitely, trading a residual lag/noise
      floor (set by s and sigma) for the ability to follow the drift. nu_op
      settles to a finite floor instead of climbing.

    Gradient of log P_{i,.}(x_{t+1} | theta) w.r.t. theta_{ij} when x_t = i:
        grad_{ij} = 1[x_{t+1} = j] - P_{ij}(theta).
    For rows other than i the gradient is 0 (observation provides no info on them
    at this step). Diagonal entry theta_{ii} held fixed = 0.

    Returns hat_theta and hat_P via softmax_rows_full(hat_theta).
    The learner accumulates all observations (no forgetting): after t updates it
    has seen N_obs = t transitions (identical for both schedules).
    """

    def __init__(self, k: int, c: float = 1.0, theta0: float = 0.0,
                 const_step: float | None = None):
        self.k = k
        self.c = c
        self.const_step = const_step
        self.theta = np.full((k, k), theta0, dtype=np.float64)
        np.fill_diagonal(self.theta, 0.0)
        self._step = 0

    def update(self, x_prev: int, x_next: int) -> None:
        self._step += 1
        if self.const_step is not None:
            eta = self.const_step
        else:
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

    def n_obs(self) -> int:
        """Number of transition observations ingested so far (= t, no window)."""
        return self._step

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
    learn_c: float = 1.0,
    const_step: float = 0.3,
) -> dict:
    """Run a single OU-drift simulation with BOTH learners on ONE trajectory.

    Both estimators (Robbins-Monro decaying step c/sqrt(t) and constant step
    `const_step`) consume the SAME realised chain x_t and see the SAME true
    transition matrix P_true(t). Only their step-size schedule differs, so the
    contrast (RM freezes vs constant step tracks) is not confounded by sampling.

    Returns a dict with arrays sampled every `measure_every` steps. Shared keys:
      times, I_opt, I_mem, n_obs.
    Per-learner keys (suffix _rm = Robbins-Monro, _const = constant step):
      I_pred_{rm,const}, eta_{rm,const}, nu_still_{rm,const}, nu_op_{rm,const},
      theta_err_{rm,const}.

    Re-equipped cost model (all in nats), identical for both learners:
      N_obs(t)    = t                         (online learner, no sliding window).
      I_mem(t)    = (K/2) ln(max(N_obs, e)).
      eta(t)      = I_pred(t) / I_mem(t)       (Landauer efficiency).
      nu_Still(t) = 1 - eta(t)                 (memory ballast, secondary).
      nu_op(t)    = 1 - I_pred(t) / I_opt(t)   (operational nostalgia, CANONICAL).
    """
    # initialise true OU logits and the two estimators (same gradient, diff step)
    theta_star = np.zeros((k, k))
    theta = init_logits(k, theta0_std, rng)
    est_rm = OnlineMLLogitEstimator(k, c=learn_c, const_step=None)
    est_const = OnlineMLLogitEstimator(k, c=learn_c, const_step=const_step)

    # initial state: stationary of current P
    P_true = softmax_rows_full(theta)
    pi0 = stationary_distribution(P_true)
    x = int(rng.choice(k, p=pi0))

    n_samples = T // measure_every
    times = np.zeros(n_samples)
    I_opt_arr = np.zeros(n_samples)
    I_mem_arr = np.zeros(n_samples)
    n_obs_arr = np.zeros(n_samples)

    # per-learner arrays
    keys_learner = ["I_pred", "eta", "nu_still", "nu_op", "theta_err"]
    arr = {f"{key}_{tag}": np.zeros(n_samples)
           for key in keys_learner for tag in ("rm", "const")}

    def measure(estimator, P_true, I_opt, I_mem):
        P_hat = estimator.estimate_P()
        I_obs = observer_predictive_info(P_true, P_hat, tau_pred)
        # finite-sample guard (NOT a DPI theorem): cap estimate at stored info; see main.ru.md §3.1.
        I_obs = min(I_obs, I_mem)
        eta = I_obs / I_mem if I_mem > 0 else 0.0
        nu_still = 1.0 - eta
        nu_op = 1.0 - (I_obs / I_opt) if I_opt > 1e-9 else 0.0
        theta_hat = estimator.estimate_theta()
        theta_err = float(np.linalg.norm(theta - theta_hat))
        return I_obs, eta, nu_still, nu_op, theta_err

    sample_idx = 0
    for t in range(1, T + 1):
        # advance OU logits
        theta = ou_step(theta, theta_star, lam, sigma, rng, dt=1.0)
        P_true = softmax_rows_full(theta)

        # advance Markov chain ONCE; both learners see the same transition
        x_next = int(rng.choice(k, p=P_true[x]))
        est_rm.update(x, x_next)
        est_const.update(x, x_next)
        x = x_next

        # measurement
        if t % measure_every == 0 and sample_idx < n_samples:
            I_opt = optimal_predictive_info(P_true, tau_pred)
            n_obs = est_rm.n_obs()  # = t (identical for both learners)
            I_mem = memory_information(n_obs, k)

            ip_rm, et_rm, ns_rm, no_rm, te_rm = measure(est_rm, P_true, I_opt, I_mem)
            ip_c, et_c, ns_c, no_c, te_c = measure(est_const, P_true, I_opt, I_mem)

            times[sample_idx] = t
            I_opt_arr[sample_idx] = I_opt
            I_mem_arr[sample_idx] = I_mem
            n_obs_arr[sample_idx] = n_obs

            arr["I_pred_rm"][sample_idx] = ip_rm
            arr["eta_rm"][sample_idx] = et_rm
            arr["nu_still_rm"][sample_idx] = ns_rm
            arr["nu_op_rm"][sample_idx] = no_rm
            arr["theta_err_rm"][sample_idx] = te_rm

            arr["I_pred_const"][sample_idx] = ip_c
            arr["eta_const"][sample_idx] = et_c
            arr["nu_still_const"][sample_idx] = ns_c
            arr["nu_op_const"][sample_idx] = no_c
            arr["theta_err_const"][sample_idx] = te_c

            sample_idx += 1

    out = dict(
        times=times[:sample_idx],
        I_opt=I_opt_arr[:sample_idx],
        I_mem=I_mem_arr[:sample_idx],
        n_obs=n_obs_arr[:sample_idx],
    )
    for key, a in arr.items():
        out[key] = a[:sample_idx]
    return out


def average_runs(runs: list[dict]) -> dict:
    shared = ["I_opt", "I_mem", "n_obs"]
    keys_learner = ["I_pred", "eta", "nu_still", "nu_op", "theta_err"]
    per_learner = [f"{key}_{tag}" for key in keys_learner for tag in ("rm", "const")]
    times = runs[0]["times"]
    out = {"times": times}
    for key in shared + per_learner:
        stacked = np.stack([r[key] for r in runs], axis=0)
        out[key] = stacked.mean(axis=0)
        out[key + "_std"] = stacked.std(axis=0)
    return out
