"""iter-006 Step 1: critical experiment.

Question: can the asymptotic of conjecture S8.1 / (S8.2) Supplementary be rescued by redefining
I_pred via a BNT-style I_infinity (cumulative Bayesian redundancy) instead of
the one-step MI used in iter-005?

Setup:
  - Same OU-drift Markov chain as in ../markov_drift_ou (k=8, K=k(k-1)=56,
    lambda=1e-3, sigma=0.01 -- SMALLER than iter-005's 0.1 to keep T<<1/sigma^2
    over the whole simulation window).
  - Two estimators of "predictive information":
      Metric A (growing-window MI):   T_w * one-step MI (closed form, Markov)
      Metric B (BNT redundancy):      L_excess(t) = sum_{s<=t} [log P_s - log hat P_s]
    See model.py for the rationale.

Main test (Metric B, BNT slope):
    Linear fit  L_excess(t)  vs  ln(lambda t)  in t in [1e3, 1e4].
    Theory: slope = K/2 = 28      (R does not enter, no budget normalisation)

Operational test (Metric B-eta):
    Linear fit  eta_L^inf(t) * t  vs  ln(lambda t)  in t in [1e3, 1e4].
    Theory: slope = K/(2R) = 28   (analogue of iter-005's failed test)

Outputs (only inside this directory; we do NOT touch paper/figs):
  - fig_iinf_growing.png            : Metric A control plot
  - fig_iinf_etat_vs_lnlambdat.png  : Metric B main fit
  - fig_iinf_Lexcess.png            : Metric B-strict (BNT redundancy) plot
  - results_summary.{txt,json}
  - run.log
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import average_runs_iinf, simulate_run_iinf

HERE = Path(__file__).resolve().parent
SUMMARY_PATH = HERE / "results_summary.txt"
SUMMARY_JSON = HERE / "results_summary.json"
LOG_PATH = HERE / "run.log"

# ----- Global parameters -----
K_STATES = 8                # number of states
K_PARAMS = K_STATES * (K_STATES - 1)   # = 56 off-diagonal logits
LAM = 1e-3                  # OU mean-reversion rate
SIGMA = 0.01                # OU diffusion amplitude (iter-006: SMALLER than iter-005)
T_TOTAL = 10_000            # simulation length in steps (= 1/sigma^2, transient bound)
N_RUNS = 30                 # independent OU realisations (start small per spec)
THETA0_STD = 0.1            # std of initial logits
MEASURE_EVERY = 100         # subsampling for measurement
R = 1.0                     # stationary energy rate
TAU_HALF = 100.0            # storage half-life (cost coefficient) -- unused in fits
LEARN_C = 3.0               # initial learning rate of online ML (Robbins-Monro c)
SEED = 20260524

# Fit window: t in [1e3, 1e4]  (per spec)
T_MIN_FIT = 1_000.0
T_MAX_FIT = 10_000.0


def linear_fit(x: np.ndarray, y: np.ndarray) -> dict:
    if x.size < 5:
        return {"ok": False, "reason": "insufficient points", "n_points": int(x.size)}
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = coef
    yhat = A @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "ok": True,
        "slope": float(slope),
        "intercept": float(intercept),
        "R2": float(r2),
        "n_points": int(x.size),
    }


def fit_in_window(times: np.ndarray, y_full: np.ndarray, lam: float,
                  t_min: float, t_max: float, positive_only: bool = False) -> dict:
    mask = (times >= t_min) & (times <= t_max)
    if positive_only:
        mask = mask & (y_full > 0)
    t = times[mask]
    if t.size < 5:
        return {"ok": False, "reason": "insufficient points", "n_points": int(t.size),
                "t_min": float(t_min), "t_max": float(t_max)}
    y = y_full[mask]
    x = np.log(np.maximum(lam * t, 1e-12))
    fit = linear_fit(x, y)
    fit["t_min"] = float(t_min)
    fit["t_max"] = float(t_max)
    return fit


def joint_fit_log_plus_linear(times: np.ndarray, y_full: np.ndarray, lam: float,
                              t_min: float, t_max: float) -> dict:
    """Fit y(t) = A * t + B * ln(lambda t) + C  on [t_min, t_max].

    Returns A (linear drift coefficient), B (BNT log coefficient), C (intercept),
    R^2 of the joint fit. Theory: B = K/2 in the BNT transient regime.
    """
    mask = (times >= t_min) & (times <= t_max)
    t = times[mask]
    if t.size < 5:
        return {"ok": False, "reason": "insufficient points", "n_points": int(t.size)}
    y = y_full[mask]
    ln = np.log(np.maximum(lam * t, 1e-12))
    A_mat = np.vstack([t, ln, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A_mat, y, rcond=None)
    A, B, C = coef
    yhat = A_mat @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "ok": True,
        "A_linear_per_t": float(A),
        "B_log_coef": float(B),
        "C_intercept": float(C),
        "R2": float(r2),
        "n_points": int(t.size),
        "t_min": float(t_min), "t_max": float(t_max),
    }


def main() -> None:
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_lines.append(msg)

    log("== iter-006 Step 1: OU + I_infinity (BNT redundancy) ==")
    log(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}")
    log(f"lambda = {LAM}, sigma = {SIGMA}, T = {T_TOTAL}, n_runs = {N_RUNS}")
    log(f"theta0_std = {THETA0_STD}, measure_every = {MEASURE_EVERY}")
    log(f"R = {R}, learn_c = {LEARN_C}")
    log(f"OU stationary variance sigma^2/(2 lambda) = {SIGMA**2 / (2 * LAM):.3e}")
    log(f"Transient-window upper bound 1/sigma^2 = {1.0 / SIGMA**2:.2e} steps")
    log(f"Fit window [{int(T_MIN_FIT)}, {int(T_MAX_FIT)}] sits entirely below transient bound: "
        f"{T_MAX_FIT <= 1.0 / SIGMA**2}")
    log(f"Theory slope (Metric B BNT-strict): K/2 = {K_PARAMS / 2.0:.1f}")
    log(f"Theory slope (Metric B operational): K/(2R) = {K_PARAMS / (2.0 * R):.1f}")

    t0 = time.time()
    runs = []
    for r in range(N_RUNS):
        rng = np.random.default_rng(SEED + r)
        run = simulate_run_iinf(
            k=K_STATES,
            T=T_TOTAL,
            lam=LAM,
            sigma=SIGMA,
            rng=rng,
            measure_every=MEASURE_EVERY,
            theta0_std=THETA0_STD,
            R=R,
            tau_half=TAU_HALF,
            learn_c=LEARN_C,
        )
        runs.append(run)
        if (r + 1) % 5 == 0:
            log(f"  done run {r + 1}/{N_RUNS} (elapsed {time.time() - t0:.1f}s)")

    log(f"Total simulation time: {time.time() - t0:.1f}s")

    avg = average_runs_iinf(runs)
    times = avg["times"]
    L_excess = avg["L_excess_cum"]
    G_obs = avg["G_obs_cum"]
    cum_I_opt = avg["cumulative_I_opt"]
    eta_inf = avg["eta_L_R_inf"]
    eta_1step = avg["eta_L_R_1step"]
    nu_arr = avg["nu"]

    # ---------- Fit 1: BNT-strict   L_excess(t)   vs   ln(lambda t) ----------
    fit_Lexcess = fit_in_window(times, L_excess, LAM, T_MIN_FIT, T_MAX_FIT)
    slope_theory_Lexcess = K_PARAMS / 2.0  # = 28
    log("\n== Quantitative analysis: BNT-strict (Metric B-strict, log-only fit) ==")
    log(f"Fit  L_excess(t)  vs  ln(lambda t)  in [{int(T_MIN_FIT)},{int(T_MAX_FIT)}]")
    log(f"  slope        = {fit_Lexcess.get('slope', float('nan')):.3f}")
    log(f"  intercept    = {fit_Lexcess.get('intercept', float('nan')):.3f}")
    log(f"  R^2          = {fit_Lexcess.get('R2', float('nan')):.3f}")
    log(f"  n_points     = {fit_Lexcess.get('n_points', 0)}")
    log(f"  theory K/2   = {slope_theory_Lexcess:.1f}")
    if fit_Lexcess.get("ok"):
        log(f"  ratio empirical/theory = {fit_Lexcess['slope'] / slope_theory_Lexcess:.3f}")

    # Joint fit:  L_excess(t) = A*t + B*ln(lam t) + C
    # In the non-stationary OU regime cumulative KL grows linearly in t (rate ~ K*sigma^2/lam)
    # plus the residual BNT redundancy (K/2) ln(lam t). Isolating the log
    # coefficient B is the cleaner BNT test.
    fit_joint = joint_fit_log_plus_linear(times, L_excess, LAM, T_MIN_FIT, T_MAX_FIT)
    log("\n== JOINT FIT: L_excess(t) = A*t + B*ln(lambda t) + C ==")
    log(f"  A (linear) = {fit_joint.get('A_linear_per_t', float('nan')):.5e}  (KL/step drift loss)")
    log(f"  B (log)    = {fit_joint.get('B_log_coef', float('nan')):.3f}")
    log(f"  C          = {fit_joint.get('C_intercept', float('nan')):.3f}")
    log(f"  R^2        = {fit_joint.get('R2', float('nan')):.4f}")
    log(f"  theory B = K/2 = {slope_theory_Lexcess:.1f}")
    if fit_joint.get("ok"):
        log(f"  ratio B/(K/2) = {fit_joint['B_log_coef'] / slope_theory_Lexcess:.3f}")

    # ---------- Fit 2: operational   eta_inf * t   vs   ln(lambda t) ----------
    # eta_L^inf(t) = G_obs(t) / (R t), so eta_inf*t = G_obs/R.
    fit_eta_inf = fit_in_window(times, eta_inf * times, LAM, T_MIN_FIT, T_MAX_FIT)
    slope_theory_eta = K_PARAMS / (2.0 * R)  # = 28
    log("\n== Quantitative analysis: operational (Metric B-eta) ==")
    log(f"Fit  eta_L^inf(t)*t  vs  ln(lambda t)  in [{int(T_MIN_FIT)},{int(T_MAX_FIT)}]")
    log(f"  slope        = {fit_eta_inf.get('slope', float('nan')):.3f}")
    log(f"  intercept    = {fit_eta_inf.get('intercept', float('nan')):.3f}")
    log(f"  R^2          = {fit_eta_inf.get('R2', float('nan')):.3f}")
    log(f"  n_points     = {fit_eta_inf.get('n_points', 0)}")
    log(f"  theory K/(2R)= {slope_theory_eta:.1f}")
    if fit_eta_inf.get("ok"):
        log(f"  ratio empirical/theory = {fit_eta_inf['slope'] / slope_theory_eta:.3f}")

    # ---------- Fit 3: control   one-step eta * t   vs   ln(lambda t) ----------
    fit_eta_1 = fit_in_window(times, eta_1step * times, LAM, T_MIN_FIT, T_MAX_FIT,
                              positive_only=True)
    log("\n== Control: one-step eta_L*t (iter-005 metric) ==")
    log(f"  slope        = {fit_eta_1.get('slope', float('nan')):.3f}")
    log(f"  R^2          = {fit_eta_1.get('R2', float('nan')):.3f}")
    if fit_eta_1.get("ok"):
        log(f"  ratio vs K/(2R)={slope_theory_eta:.0f}: "
            f"{fit_eta_1['slope'] / slope_theory_eta:.3f}")

    # ---------- Diagnostics ----------
    late_mask = times > 0.5 * T_TOTAL
    nu_late = float(np.mean(nu_arr[late_mask]))
    nu_late_std = float(np.std(nu_arr[late_mask]))
    I_opt_late = float(np.mean(avg["I_opt_1step"][late_mask]))
    log("\n== Diagnostics ==")
    log(f"Late-time nu = {nu_late:.3f} +- {nu_late_std:.3f}")
    log(f"Late-time <I_opt_1step> = {I_opt_late:.3e}")
    log(f"<theta_err>_late = {float(np.mean(avg['theta_err'][late_mask])):.3f}")
    log(f"L_excess at t={int(times[-1])}: {L_excess[-1]:.3f}")
    log(f"G_obs    at t={int(times[-1])}: {G_obs[-1]:.3f}")
    log(f"sum_s I_opt at t={int(times[-1])}: {cum_I_opt[-1]:.3f}")

    # ---------- FIG A: growing-window metric ----------
    growing_T = avg["growing_T"]
    growing_I_obs = avg["growing_I_obs_T"]
    growing_I_obs_std = avg["growing_I_obs_T_std"]
    growing_I_opt = avg["growing_I_opt_T"]
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.errorbar(growing_T, growing_I_obs, yerr=growing_I_obs_std,
                fmt="o-", color="tab:blue", lw=1.5, ms=5,
                label=r"Metric A: $I_{obs}^{A}(T_f, T_w) = T_w \cdot I_1^{obs}$")
    ax.plot(growing_T, growing_I_opt, "s--", color="tab:orange", lw=1.2, ms=4,
            label=r"oracle: $T_w \cdot I_1^{opt}$")
    # theory line for I_infinity ~ (K/2) ln(lambda T_w)  (BNT 2001 Class II)
    bnt_curve = (K_PARAMS / 2.0) * np.log(np.maximum(LAM * growing_T, 1e-12))
    bnt_curve_shifted = bnt_curve - bnt_curve.min() + growing_I_opt.min()  # for visual
    ax.plot(growing_T, bnt_curve_shifted, ":", color="tab:red", lw=1.2,
            label=r"BNT (S8.2) target $(K/2)\ln(\lambda T_w)$ (anchored)")
    ax.set_xscale("log")
    ax.set_xlabel(r"window size $T_w$ (steps each side of $t=T_{\mathrm{final}}$)")
    ax.set_ylabel("predictive information (nats)")
    ax.set_title(rf"Metric A (growing-window MI) at $t = T_{{\rm final}} = {T_TOTAL}$, "
                 rf"$N_{{\rm runs}}={N_RUNS}$")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    figA_path = HERE / "fig_iinf_growing.png"
    fig.savefig(figA_path, dpi=150)
    plt.close(fig)
    log(f"saved {figA_path}")

    # ---------- FIG B-eta: main operational fit ----------
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    x_all = np.log(np.maximum(LAM * times, 1e-12))
    y_all = eta_inf * times
    ax.plot(x_all, y_all, "o", ms=3.5, alpha=0.6, color="tab:blue",
            label=r"data: $\eta_L^\infty(t)\cdot t = G_{\rm obs}(t)/R$")
    mask_fit = (times >= T_MIN_FIT) & (times <= T_MAX_FIT)
    if fit_eta_inf.get("ok"):
        x_fit = np.linspace(x_all[mask_fit].min(), x_all[mask_fit].max(), 200)
        y_fit = fit_eta_inf["slope"] * x_fit + fit_eta_inf["intercept"]
        ax.plot(x_fit, y_fit, "-", color="tab:blue", lw=1.6,
                label=rf"linear fit: slope = {fit_eta_inf['slope']:.2f} "
                      rf"($R^2$={fit_eta_inf['R2']:.2f})")
        # theory anchor
        x_mid = 0.5 * (x_all[mask_fit].min() + x_all[mask_fit].max())
        y_mid = float(np.mean(y_all[mask_fit]))
        b_th = y_mid - slope_theory_eta * x_mid
        ax.plot(x_fit, slope_theory_eta * x_fit + b_th, "--", color="tab:red", lw=1.4,
                label=rf"theory slope $K/(2R)={slope_theory_eta:.0f}$ (anchored)")
    ax.axvspan(np.log(LAM * T_MIN_FIT), np.log(LAM * T_MAX_FIT),
               color="grey", alpha=0.10, label=f"fit window [{int(T_MIN_FIT)},{int(T_MAX_FIT)}]")
    ax.set_xlabel(r"$\ln(\lambda t)$")
    ax.set_ylabel(r"$\eta_L^\infty(t)\cdot t\ \ \mathrm{[nats]}$")
    ax.set_title(r"Main fit (Metric B operational): $\eta_L^\infty\!\cdot\!t$ vs $\ln(\lambda t)$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    figBe_path = HERE / "fig_iinf_etat_vs_lnlambdat.png"
    fig.savefig(figBe_path, dpi=150)
    plt.close(fig)
    log(f"saved {figBe_path}")

    # ---------- FIG B-strict: BNT redundancy plot ----------
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(x_all, L_excess, "o", ms=3.5, alpha=0.6, color="tab:green",
            label=r"data: $L_{\rm excess}(t)$  (cumulative BNT redundancy)")
    if fit_Lexcess.get("ok"):
        x_fit = np.linspace(x_all[mask_fit].min(), x_all[mask_fit].max(), 200)
        y_fit = fit_Lexcess["slope"] * x_fit + fit_Lexcess["intercept"]
        ax.plot(x_fit, y_fit, "-", color="tab:green", lw=1.6,
                label=rf"linear fit: slope = {fit_Lexcess['slope']:.2f} "
                      rf"($R^2$={fit_Lexcess['R2']:.2f})")
        x_mid = 0.5 * (x_all[mask_fit].min() + x_all[mask_fit].max())
        y_mid = float(np.mean(L_excess[mask_fit]))
        b_th = y_mid - slope_theory_Lexcess * x_mid
        ax.plot(x_fit, slope_theory_Lexcess * x_fit + b_th, "--", color="tab:red", lw=1.4,
                label=rf"BNT theory slope $K/2={slope_theory_Lexcess:.0f}$ (anchored)")
    ax.axvspan(np.log(LAM * T_MIN_FIT), np.log(LAM * T_MAX_FIT),
               color="grey", alpha=0.10, label=f"fit window [{int(T_MIN_FIT)},{int(T_MAX_FIT)}]")
    ax.set_xlabel(r"$\ln(\lambda t)$")
    ax.set_ylabel(r"$L_{\rm excess}(t)\ \ \mathrm{[nats]}$")
    ax.set_title("BNT-strict redundancy: $L_{\\rm excess}(t)$ vs $\\ln(\\lambda t)$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    figBs_path = HERE / "fig_iinf_Lexcess.png"
    fig.savefig(figBs_path, dpi=150)
    plt.close(fig)
    log(f"saved {figBs_path}")

    # ---------- Write summary ----------
    summary = {
        "experiment": "iter-006 step 1: OU + I_infinity (BNT redundancy)",
        "k": K_STATES,
        "K_params": K_PARAMS,
        "lambda": LAM,
        "sigma": SIGMA,
        "T_total": T_TOTAL,
        "n_runs": N_RUNS,
        "theta0_std": THETA0_STD,
        "R": R,
        "tau_half": TAU_HALF,
        "learn_c": LEARN_C,
        "transient_window_upper": 1.0 / SIGMA**2,
        "fit_window": [T_MIN_FIT, T_MAX_FIT],
        "slope_theory_K_over_2": K_PARAMS / 2.0,
        "slope_theory_K_over_2R": K_PARAMS / (2.0 * R),
        "fit_Lexcess_vs_log_lam_t": fit_Lexcess,
        "fit_Lexcess_joint_linear_plus_log": fit_joint,
        "fit_eta_inf_t_vs_log_lam_t": fit_eta_inf,
        "fit_eta_1step_t_vs_log_lam_t_control": fit_eta_1,
        "late_time": {
            "nu_mean": nu_late, "nu_std": nu_late_std,
            "I_opt_1step_mean": I_opt_late,
            "theta_err_mean": float(np.mean(avg["theta_err"][late_mask])),
        },
        "final": {
            "t": int(times[-1]),
            "L_excess": float(L_excess[-1]),
            "G_obs": float(G_obs[-1]),
            "cum_I_opt": float(cum_I_opt[-1]),
        },
        "metric_A_growing": {
            "T_w": [int(v) for v in avg["growing_T"]],
            "I_obs": [float(v) for v in avg["growing_I_obs_T"]],
            "I_opt": [float(v) for v in avg["growing_I_opt_T"]],
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    with SUMMARY_PATH.open("w", encoding="utf-8") as f:
        f.write("=== iter-006 Step 1: OU + I_infinity (BNT redundancy) ===\n")
        f.write(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}\n")
        f.write(f"lambda = {LAM}, sigma = {SIGMA}, T = {T_TOTAL}, n_runs = {N_RUNS}\n")
        f.write(f"OU stat. variance sigma^2/(2 lambda) = {SIGMA**2 / (2 * LAM):.3e}\n")
        f.write(f"Transient bound 1/sigma^2 = {1.0 / SIGMA**2:.2e}\n")
        f.write(f"Fit window [{int(T_MIN_FIT)}, {int(T_MAX_FIT)}]  -- within transient bound: "
                f"{T_MAX_FIT <= 1.0 / SIGMA**2}\n\n")
        f.write("--- Metric B-strict (log-only): L_excess(t) vs ln(lambda t) ---\n")
        f.write(f"  slope     = {fit_Lexcess.get('slope', float('nan')):.3f}\n")
        f.write(f"  intercept = {fit_Lexcess.get('intercept', float('nan')):.3f}\n")
        f.write(f"  R^2       = {fit_Lexcess.get('R2', float('nan')):.3f}\n")
        f.write(f"  theory K/2 = {K_PARAMS / 2.0:.1f}\n")
        if fit_Lexcess.get("ok"):
            f.write(f"  ratio empirical/theory = {fit_Lexcess['slope'] / (K_PARAMS / 2.0):.3f}\n")
        f.write("\n--- Joint fit: L_excess(t) = A*t + B*ln(lambda t) + C ---\n")
        f.write(f"  A (linear drift) = {fit_joint.get('A_linear_per_t', float('nan')):.5e}\n")
        f.write(f"  B (log coef)     = {fit_joint.get('B_log_coef', float('nan')):.3f}\n")
        f.write(f"  C                = {fit_joint.get('C_intercept', float('nan')):.3f}\n")
        f.write(f"  R^2              = {fit_joint.get('R2', float('nan')):.4f}\n")
        f.write(f"  theory B = K/2   = {K_PARAMS / 2.0:.1f}\n")
        if fit_joint.get("ok"):
            f.write(f"  ratio B/(K/2)   = {fit_joint['B_log_coef'] / (K_PARAMS / 2.0):.3f}\n")
        f.write("\n--- Metric B-eta (operational): eta_L^inf(t)*t vs ln(lambda t) ---\n")
        f.write(f"  slope     = {fit_eta_inf.get('slope', float('nan')):.3f}\n")
        f.write(f"  R^2       = {fit_eta_inf.get('R2', float('nan')):.3f}\n")
        f.write(f"  theory K/(2R) = {K_PARAMS / (2.0 * R):.1f}\n")
        if fit_eta_inf.get("ok"):
            f.write(f"  ratio empirical/theory = {fit_eta_inf['slope'] / (K_PARAMS / (2.0 * R)):.3f}\n")
        f.write("\n--- Control: one-step eta*t (iter-005 metric) ---\n")
        f.write(f"  slope     = {fit_eta_1.get('slope', float('nan')):.3f}\n")
        f.write(f"  R^2       = {fit_eta_1.get('R2', float('nan')):.3f}\n")
        f.write("\nLate-time diagnostics:\n")
        f.write(f"  nu = {nu_late:.3f} +- {nu_late_std:.3f}\n")
        f.write(f"  <I_opt_1step> = {I_opt_late:.3e}\n")
        f.write(f"  L_excess({int(times[-1])}) = {L_excess[-1]:.3f}\n")
        f.write(f"  G_obs({int(times[-1])})    = {G_obs[-1]:.3f}\n")

    LOG_PATH.write_text("\n".join(log_lines), encoding="utf-8")
    print(f"\nWrote {SUMMARY_PATH}, {SUMMARY_JSON}, {LOG_PATH}")


if __name__ == "__main__":
    main()
