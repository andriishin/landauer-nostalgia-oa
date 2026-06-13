"""Model for iter-006 Step 1: OU-drift Markov chain with I_infinity estimators.

This file re-uses the OU/softmax/learner core from ``../markov_drift_ou/model.py``
and adds the two excess-entropy / predictive-information estimators required to
test whether replacing the one-step MI in conjecture S8.1 / (S8.2) Supplementary
with a BNT-style I_infinity
recovers the predicted (K/2) ln(lambda t) growth.

Two estimators of "predictive information" are provided:

Metric A (growing-window MI):
    I_pred^A(t, T) := I(B_t^- ; B_t^+),
    where B_t^- = (X_{t-T+1},...,X_t) and B_t^+ = (X_{t+1},...,X_{t+T}) are
    half-blocks of length T sitting either side of time t. We do NOT estimate the
    high-dimensional joint distribution over 8^T symbols (intractable for T>~3);
    instead, because the chain is (slowly) Markovian we use the analytical
    decomposition for a known transition matrix:

        I(B^- ; B^+) = T * I_1 + O(1),

    where I_1 = I(X_t ; X_{t+1}) is the one-step MI under the current P(t). For
    a *learner* with estimate hat_P(t), the observable analogue is

        I_obs^A(t, T) = T * I_obs(P(t), hat_P(t), tau=1),

    which is exactly T times the one-step observer MI used in iter-005. This
    rescaling DOES NOT cure the structural ceiling, because for a fixed truth P
    the per-symbol mutual information is bounded by ln k, and the time-averaged
    quantity is therefore still O(ln k). We compute and report it as a control
    that makes the ceiling explicit.

Metric B (Bayesian information gain / accumulated log-likelihood):
    G_obs(t) := sum_{s<=t} [log hat_P(X_s | X_{s-1}) - log(1/k)]
              = sum_{s<=t} [log hat_P(X_s | X_{s-1}) + log k].

    This is the cumulative log-likelihood improvement of the learner's model
    over a uniform-output baseline. By the chain rule for KL,

        E G_obs(t) = sum_{s<=t} [ H_max - H(X_s|X_{s-1}) - KL(P_s || hat_P_s) ]
                   = sum_{s<=t} I_opt(P_s, 1) - cumulative_KL_loss(t).

    For an ideal learner the KL loss is O((K/2) ln t) (the standard Bayesian
    "predictive description length" theorem; see BNT 2001 eq. for Class II
    regular models), so G_obs grows linearly in t with a *logarithmic correction*
    of size (K/2) ln(lambda t). After dividing by N_max(t) = R*t we obtain

        eta_L^infinity(t) := G_obs(t) / (R t)  -> I_opt_inf  - (K/(2R)) ln(lam t)/t,

    so the quantity  (I_opt_inf * t - G_obs(t))  vs  ln(lam t) should have
    SLOPE = K/(2R) = 28 in the transient regime.

This file also provides a *cleaner* version of the same idea that doesn't need
I_opt_inf: the cumulative excess code-length

        L_excess(t) := sum_{s<=t} [ -log hat_P(X_s|X_{s-1}) + log P(X_s|X_{s-1}) ]
                     = sum_{s<=t} [ log P_s(X_s|X_{s-1}) - log hat_P_s(X_s|X_{s-1}) ],

    which is the *one-sample* analogue of cumulative KL(P || hat_P). For an
    ideal Bayesian learner with K-dimensional regular parameter, BNT 2001
    predicts  E L_excess(t) ~ (K/2) ln(lambda t) + O(1)  (this is the universal
    Bayesian redundancy, transported into the slowly-drifting OU regime by the
    transient bound t << 1/sigma^2).

    The "BNT slope test" is then simply:

        slope of  L_excess(t)  vs  ln(lambda t)   ?=   K/2 = 28   (R does not enter)

    No budget normalisation is needed; this is the cleanest possible test of
    "does the one-step MI ceiling go away when we measure cumulative redundancy
    instead". We compute and plot this as well.

Sanity: the simulation runs at lambda = 1e-3, sigma = 0.01 (much smaller than
the iter-005 sigma=0.1) so that the transient window 1/sigma^2 = 1e4 is the
*entire simulation length*; this is the regime BNT eq. for slowly-drifting
parameters formally addresses.

All entropies / KL in nats. Reuses helpers from the iter-005 OU module.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

# Load the iter-005 OU module directly from its file path under a unique module
# name so we re-use its OU step, softmax row normaliser, KL helpers, and online
# learner without collision with our own local `model` module.
_OU_PARENT = Path(__file__).resolve().parent.parent / "markov_drift_ou"
_OU_MODULE_PATH = _OU_PARENT / "model.py"
_spec = importlib.util.spec_from_file_location("markov_drift_ou_model", _OU_MODULE_PATH)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Cannot load iter-005 OU model from {_OU_MODULE_PATH}")
_ou = importlib.util.module_from_spec(_spec)
sys.modules["markov_drift_ou_model"] = _ou
_spec.loader.exec_module(_ou)

OnlineMLLogitEstimator = _ou.OnlineMLLogitEstimator
entropy_of_distribution = _ou.entropy_of_distribution
init_logits = _ou.init_logits
kl_div_rows = _ou.kl_div_rows
matrix_power = _ou.matrix_power
optimal_predictive_info = _ou.optimal_predictive_info
observer_predictive_info = _ou.observer_predictive_info
ou_step = _ou.ou_step
softmax_rows_full = _ou.softmax_rows_full
stationary_distribution = _ou.stationary_distribution


def per_symbol_log_density(P: np.ndarray, x_prev: int, x_next: int,
                            eps: float = 1e-12) -> float:
    """log P(x_next | x_prev) under transition matrix P."""
    return float(np.log(max(P[x_prev, x_next], eps)))


def simulate_run_iinf(
    k: int,
    T: int,
    lam: float,
    sigma: float,
    rng: np.random.Generator,
    measure_every: int = 100,
    theta0_std: float = 0.1,
    R: float = 1.0,
    N0: float = 0.0,
    tau_half: float = 100.0,
    learn_c: float = 1.0,
    growing_T_grid: tuple[int, ...] = (10, 30, 100, 300, 1000, 3000),
) -> dict:
    """One OU realisation, recording I_infinity-style functionals.

    Returns dict containing, sampled every ``measure_every`` steps:
      times             : t at each measurement
      I_pred_1step      : observer one-step MI at this t (iter-005 metric, for control)
      I_opt_1step       : optimal one-step MI under true P(t)
      L_excess_cum      : cumulative excess code-length sum_{s<=t} [log P_s - log hat P_s]
                          evaluated on realised samples (Metric B-strict, BNT redundancy)
      G_obs_cum         : cumulative info gain  sum_{s<=t} [log hat P_s + log k]
                          (Metric B "operational" form, vs uniform baseline)
      cumulative_I_opt  : sum_{s<=t} I_opt_1step(s) (oracle counterpart of G_obs)
      eta_L_R_1step     : I_pred_1step / (R t)  (iter-005 efficiency, control)
      eta_L_inf_R       : G_obs_cum  / (R t)    (Metric B operational efficiency)
      nu                : 1 - I_pred_1step / I_opt_1step  (instantaneous)
      theta_err         : ||theta(t) - hat_theta(t)||_F (posterior concentration diag.)

    And, only at the final time step (T), the growing-window metric A:
      growing_T            : array of window sizes T_w
      growing_I_obs_T      : I_obs^A(T_final, T_w)  = T_w * one-step observer MI
      growing_I_opt_T      : T_w * optimal one-step MI
    (these arrays are scalars-per-T_w, repeated at every measurement time would
    be O(N_samples * N_T_grid) calls which is wasteful; we just record final t.)
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
    I_pred_1step = np.zeros(n_samples)
    I_opt_1step = np.zeros(n_samples)
    L_excess_cum_arr = np.zeros(n_samples)
    G_obs_cum_arr = np.zeros(n_samples)
    cumulative_I_opt_arr = np.zeros(n_samples)
    nu_arr = np.zeros(n_samples)
    NmaxR_arr = np.zeros(n_samples)
    eta_R_1step_arr = np.zeros(n_samples)
    eta_R_inf_arr = np.zeros(n_samples)
    theta_err_arr = np.zeros(n_samples)

    # also record the trajectory of states + the (changing) true P at each step;
    # we need them for the BNT redundancy update. We avoid storing the entire
    # sequence by accumulating online.
    L_excess_cum = 0.0
    G_obs_cum = 0.0
    cumulative_I_opt = 0.0
    log_k = float(np.log(k))

    sample_idx = 0
    for t in range(1, T + 1):
        # advance OU logits
        theta = ou_step(theta, theta_star, lam, sigma, rng, dt=1.0)
        P_true = softmax_rows_full(theta)

        # advance Markov chain
        x_next = int(rng.choice(k, p=P_true[x]))

        # PRE-UPDATE log-density under current learner: this is what the
        # learner predicted *before* seeing x_next (i.e. uses hat_theta(t-1))
        P_hat = estimator.estimate_P()
        log_p_true = per_symbol_log_density(P_true, x, x_next)
        log_p_hat = per_symbol_log_density(P_hat, x, x_next)
        L_excess_cum += (log_p_true - log_p_hat)
        G_obs_cum += (log_p_hat + log_k)

        # one-step optimal MI on the true current P (instantaneous oracle)
        I_opt_step = optimal_predictive_info(P_true, tau=1)
        cumulative_I_opt += I_opt_step

        # learner update (online)
        estimator.update(x, x_next)
        x = x_next

        # measurement
        if t % measure_every == 0 and sample_idx < n_samples:
            # snapshot one-step observer MI (post-update P_hat at this t)
            P_hat_post = estimator.estimate_P()
            I_obs_step = observer_predictive_info(P_true, P_hat_post, tau=1)
            nu = 1.0 - (I_obs_step / I_opt_step) if I_opt_step > 1e-9 else 0.0
            NmaxR = N0 + R * t
            theta_err = float(np.linalg.norm(theta - estimator.estimate_theta()))

            times[sample_idx] = t
            I_pred_1step[sample_idx] = I_obs_step
            I_opt_1step[sample_idx] = I_opt_step
            L_excess_cum_arr[sample_idx] = L_excess_cum
            G_obs_cum_arr[sample_idx] = G_obs_cum
            cumulative_I_opt_arr[sample_idx] = cumulative_I_opt
            nu_arr[sample_idx] = nu
            NmaxR_arr[sample_idx] = NmaxR
            eta_R_1step_arr[sample_idx] = I_obs_step / NmaxR if NmaxR > 0 else 0.0
            eta_R_inf_arr[sample_idx] = G_obs_cum / NmaxR if NmaxR > 0 else 0.0
            theta_err_arr[sample_idx] = theta_err
            sample_idx += 1

    # ---------------- Metric A: growing-window MI at final t ----------------
    # I(B^-_t ; B^+_t) for a Markov chain decomposes to a sum of one-step MIs.
    # For a slowly-varying P, with window length T_w on each side and current
    # truth P_T, we use the closed form  I_block = T_w * I_1 + O(1).
    # The observer counterpart with hat P_T replaces I_1 by I_obs_1.
    # This is faithful to the analytical content of "predictive information at
    # horizon T_w" and avoids estimating a 8^T_w-dimensional joint.
    P_hat_final = estimator.estimate_P()
    I_opt_1_final = optimal_predictive_info(P_true, tau=1)
    I_obs_1_final = observer_predictive_info(P_true, P_hat_final, tau=1)
    growing_T_arr = np.array(growing_T_grid, dtype=float)
    growing_I_obs_T = growing_T_arr * I_obs_1_final
    growing_I_opt_T = growing_T_arr * I_opt_1_final

    return dict(
        times=times[:sample_idx],
        I_pred_1step=I_pred_1step[:sample_idx],
        I_opt_1step=I_opt_1step[:sample_idx],
        L_excess_cum=L_excess_cum_arr[:sample_idx],
        G_obs_cum=G_obs_cum_arr[:sample_idx],
        cumulative_I_opt=cumulative_I_opt_arr[:sample_idx],
        nu=nu_arr[:sample_idx],
        N_max_R=NmaxR_arr[:sample_idx],
        eta_L_R_1step=eta_R_1step_arr[:sample_idx],
        eta_L_R_inf=eta_R_inf_arr[:sample_idx],
        theta_err=theta_err_arr[:sample_idx],
        growing_T=growing_T_arr,
        growing_I_obs_T=growing_I_obs_T,
        growing_I_opt_T=growing_I_opt_T,
    )


def average_runs_iinf(runs: list[dict]) -> dict:
    keys = [
        "I_pred_1step", "I_opt_1step", "L_excess_cum", "G_obs_cum",
        "cumulative_I_opt", "nu", "N_max_R", "eta_L_R_1step", "eta_L_R_inf",
        "theta_err",
    ]
    times = runs[0]["times"]
    out = {"times": times}
    for key in keys:
        stacked = np.stack([r[key] for r in runs], axis=0)
        out[key] = stacked.mean(axis=0)
        out[key + "_std"] = stacked.std(axis=0)

    # growing-window metric A: average over runs (scalar per T_w)
    out["growing_T"] = runs[0]["growing_T"]
    stacked_obs = np.stack([r["growing_I_obs_T"] for r in runs], axis=0)
    stacked_opt = np.stack([r["growing_I_opt_T"] for r in runs], axis=0)
    out["growing_I_obs_T"] = stacked_obs.mean(axis=0)
    out["growing_I_obs_T_std"] = stacked_obs.std(axis=0)
    out["growing_I_opt_T"] = stacked_opt.mean(axis=0)
    out["growing_I_opt_T_std"] = stacked_opt.std(axis=0)

    return out
