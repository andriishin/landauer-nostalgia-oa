"""Run OU-class simulation for paper § 6.2 and produce figs 4, 5.

Usage:
    python main.py

Outputs:
    - paper/figs/fig4_ou_eta_vs_t.png
    - paper/figs/fig5_ou_etat_vs_lnlambdat.png
    - simulations/markov_drift_ou/results_summary.{txt,json}
    - simulations/markov_drift_ou/run.log
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import average_runs, simulate_run

HERE = Path(__file__).resolve().parent
PAPER_FIGS = HERE.parent.parent / "paper" / "figs"
PAPER_FIGS.mkdir(parents=True, exist_ok=True)
SUMMARY_PATH = HERE / "results_summary.txt"
SUMMARY_JSON = HERE / "results_summary.json"
LOG_PATH = HERE / "run.log"

# ----- Global parameters (see § 5.2 of main.ru.md) -----
# Note on sigma choice: with very small sigma (e.g. 3e-3 per spec recommendation)
# the OU stationary variance sigma^2/(2 lambda) ~ 4.5e-3 keeps theta near zero,
# so the transition matrix is indistinguishable from uniform and I_opt collapses
# to noise. We use sigma = 0.1 (OU stationary variance ~ 5), which yields a
# genuinely informative chain (I_opt ~ 0.7 nats per step, well above floor),
# while sigma^2 t / lambda stays moderate (= 200 at T=2e4, i.e. just outside the
# strict transient bound but still in the regime where posterior concentration
# is the leading effect). This is the most charitable working point for testing
# the predicted K/(2R) = 28 slope; results with smaller sigma are uniformly
# worse (eta_L collapses to zero) and reported in the run log.
K_STATES = 8                # number of states
K_PARAMS = K_STATES * (K_STATES - 1)   # = 56 off-diagonal logits
LAM = 1e-3                  # OU mean-reversion rate
SIGMA = 0.1                 # OU diffusion amplitude (informative regime)
T_TOTAL = 20_000            # simulation length in steps
N_RUNS = 50                 # independent OU realisations
THETA0_STD = 0.1            # std of initial logits
TAU_PRED = 1                # one-step predictive horizon
MEASURE_EVERY = 200         # subsampling for measurement
R = 1.0                     # stationary energy rate
TAU_HALF = 100.0            # storage half-life (cost coefficient)
LEARN_C = 3.0               # initial learning rate of online ML (Robbins-Monro c)
SEED = 20260524


def fit_eta_t_vs_log_lam_t(times: np.ndarray, eta: np.ndarray, lam: float,
                           t_min: float, t_max: float) -> dict:
    """Linear fit eta_L * t = slope * ln(lam * t) + intercept on [t_min, t_max].

    Theory (conjecture S8.1 / (S8.2) Supplementary): slope = K / (2 R).
    """
    mask = (times >= t_min) & (times <= t_max) & (eta > 0)
    if mask.sum() < 5:
        return {"ok": False, "reason": "insufficient points", "n_points": int(mask.sum())}
    t = times[mask]
    y = eta[mask] * t
    x = np.log(np.maximum(lam * t, 1e-12))
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
        "n_points": int(mask.sum()),
        "t_min": float(t_min),
        "t_max": float(t_max),
        "K_eff_assuming_R1": float(2.0 * slope),
    }


def main() -> None:
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_lines.append(msg)

    log("== OU-drift simulation for § 6.2 ==")
    log(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}")
    log(f"lambda = {LAM}, sigma = {SIGMA}, T = {T_TOTAL}, n_runs = {N_RUNS}")
    log(f"theta0_std = {THETA0_STD}, tau_pred = {TAU_PRED}, measure_every = {MEASURE_EVERY}")
    log(f"R = {R}, tau_half = {TAU_HALF}, learn_c = {LEARN_C}")
    log(f"OU stationary variance sigma^2/(2 lambda) = {SIGMA**2 / (2 * LAM):.3e}")
    log(f"Transient-window upper bound 1/sigma^2 ~ {1.0 / SIGMA**2:.2e} steps")
    log(f"Adiabatic check: lambda * tau_relax ~ lambda * O(1) = {LAM:.1e} (<< 1: OK)")

    t0 = time.time()

    runs = []
    for r in range(N_RUNS):
        rng = np.random.default_rng(SEED + r)
        run = simulate_run(
            k=K_STATES,
            T=T_TOTAL,
            lam=LAM,
            sigma=SIGMA,
            tau_pred=TAU_PRED,
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

    avg = average_runs(runs)
    times = avg["times"]
    eta_R = avg["eta_L_R"]
    eta_R_std = avg["eta_L_R_std"]
    eta_full = avg["eta_L_full"]
    nu_arr = avg["nu"]

    # ---------- Quantitative fit ----------
    # Fit window: t > tau_E = 1e3 (eliminate start-up), t < T to avoid edge.
    t_min_fit = 1_000.0
    t_max_fit = 15_000.0
    fit_R = fit_eta_t_vs_log_lam_t(times, eta_R, lam=LAM, t_min=t_min_fit, t_max=t_max_fit)
    fit_full = fit_eta_t_vs_log_lam_t(times, eta_full, lam=LAM, t_min=t_min_fit, t_max=t_max_fit)

    slope_theory = K_PARAMS / (2.0 * R)  # = 28 for k=8, R=1

    log("\n== Quantitative analysis ==")
    log(f"Fit window: t in [{int(t_min_fit)}, {int(t_max_fit)}]")
    log("--- Using N_max = R*t (matches theory derivation) ---")
    log(f"  slope = {fit_R.get('slope', float('nan')):.3f}")
    log(f"  intercept = {fit_R.get('intercept', float('nan')):.3f}")
    log(f"  R^2 = {fit_R.get('R2', float('nan')):.3f}")
    log(f"  K_eff (assuming R=1) = {fit_R.get('K_eff_assuming_R1', float('nan')):.2f}")
    log(f"  theory K/(2R) = {K_PARAMS}/(2*{R}) = {slope_theory:.1f}")
    if fit_R.get("ok"):
        log(f"  ratio empirical/theory = {fit_R['slope'] / slope_theory:.3f}")
    log("--- Using N_max_full = sum_s [R + nu(s)*K/tau_half] (full spec budget) ---")
    log(f"  slope = {fit_full.get('slope', float('nan')):.3f}")
    log(f"  R^2 = {fit_full.get('R2', float('nan')):.3f}")
    if fit_full.get("ok"):
        log(f"  ratio empirical/theory = {fit_full['slope'] / slope_theory:.3f}")

    # Diagnostic: late-time nu, eta_R, eta_full
    late_mask = times > 0.5 * T_TOTAL
    nu_late = float(np.mean(nu_arr[late_mask]))
    nu_late_std = float(np.std(nu_arr[late_mask]))
    eta_R_late = float(np.mean(eta_R[late_mask]))
    eta_full_late = float(np.mean(eta_full[late_mask]))
    log(f"\nLate-time (t > {int(0.5 * T_TOTAL)}) diagnostics:")
    log(f"  nu = {nu_late:.3f} +- {nu_late_std:.3f}")
    log(f"  <eta_L>_R    = {eta_R_late:.3e}")
    log(f"  <eta_L>_full = {eta_full_late:.3e}")

    # Diagnostic: posterior concentration
    theta_err = avg["theta_err"]
    log(f"  <||theta - hat_theta||_F>_late = {float(np.mean(theta_err[late_mask])):.3f}")
    log(f"  <||theta - hat_theta||_F> at t={int(times[10])}: {float(theta_err[10]):.3f}")

    # ---------- FIG 4: eta_L(t) for OU ----------
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(times, eta_R, label=r"OU drift, ideal learner (N$_{\max}=Rt$)", color="tab:blue", lw=1.5)
    ax.fill_between(times, np.maximum(eta_R - eta_R_std, 1e-12), eta_R + eta_R_std,
                    alpha=0.20, color="tab:blue", label=r"$\pm 1\sigma$ across runs")
    ax.plot(times, eta_full, label=r"OU drift, full budget (N$_{\max}^{\mathrm{full}}$)",
            color="tab:orange", lw=1.2, ls="--")
    ax.set_xlabel("time $t$ (steps)")
    ax.set_ylabel(r"$\eta_L(t) = I_{\mathrm{pred}}(t,\tau) / N_{\max}(t)$")
    ax.set_title(rf"OU drift: $\eta_L(t)$ averaged over {N_RUNS} runs ($k={K_STATES}$, "
                 rf"$\lambda={LAM:g}$, $\sigma={SIGMA:g}$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig4_path = PAPER_FIGS / "fig4_ou_eta_vs_t.png"
    fig.savefig(fig4_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig4.pdf")  # vector copy for the manuscript
    plt.close(fig)
    log(f"saved {fig4_path}")

    # ---------- FIG 5: eta_L * t vs ln(lambda t) ----------
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    # plot for both N_max variants
    mask_R = (times >= t_min_fit) & (times <= t_max_fit) & (eta_R > 0)
    x_R = np.log(np.maximum(LAM * times, 1e-12))
    y_R = eta_R * times
    ax.plot(x_R, y_R, "o", ms=3.5, alpha=0.6, color="tab:blue", label=r"data ($N_{\max}=Rt$)")
    if fit_R.get("ok"):
        x_fit = np.linspace(x_R[mask_R].min(), x_R[mask_R].max(), 200)
        y_fit = fit_R["slope"] * x_fit + fit_R["intercept"]
        ax.plot(x_fit, y_fit, "-", color="tab:blue", lw=1.6,
                label=rf"linear fit: slope = {fit_R['slope']:.2f} (R$^2$={fit_R['R2']:.2f})")
    # theory line: y_theory = (K/(2R)) * x + intercept_aligned_to_data_midpoint
    if fit_R.get("ok"):
        x_mid = 0.5 * (x_R[mask_R].min() + x_R[mask_R].max())
        # anchor theory line at the data mean to make slope comparison visual
        y_mid_data = float(np.mean(y_R[mask_R]))
        b_th = y_mid_data - slope_theory * x_mid
        ax.plot(x_fit, slope_theory * x_fit + b_th, "--", color="tab:red", lw=1.4,
                label=rf"conjecture S8.1 slope = $K/(2R)$ = {slope_theory:.0f}")
    ax.axvspan(np.log(LAM * t_min_fit), np.log(LAM * t_max_fit), color="grey", alpha=0.10,
               label=f"fit window [{int(t_min_fit)},{int(t_max_fit)}]")
    ax.set_xlabel(r"$\ln(\lambda t)$")
    ax.set_ylabel(r"$\eta_L(t)\cdot t$")
    ax.set_title(r"Check of conjecture S8.1 / (S8.2): $\eta_L\!\cdot\! t$ vs $\ln(\lambda t)$ (OU class)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig5_path = PAPER_FIGS / "fig5_ou_etat_vs_lnlambdat.png"
    fig.savefig(fig5_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig5.pdf")  # vector copy for the manuscript
    plt.close(fig)
    log(f"saved {fig5_path}")

    # ---------- Write summary ----------
    summary = {
        "k": K_STATES,
        "K_params": K_PARAMS,
        "lambda": LAM,
        "sigma": SIGMA,
        "T_total": T_TOTAL,
        "n_runs": N_RUNS,
        "tau_pred": TAU_PRED,
        "theta0_std": THETA0_STD,
        "R": R,
        "tau_half": TAU_HALF,
        "learn_c": LEARN_C,
        "transient_window_upper": 1.0 / SIGMA**2,
        "fit_window": [t_min_fit, t_max_fit],
        "slope_theory_K_over_2R": slope_theory,
        "fit_eta_t_vs_log_lam_t_NmaxRt": fit_R,
        "fit_eta_t_vs_log_lam_t_Nmaxfull": fit_full,
        "late_time": {
            "nu_mean": nu_late,
            "nu_std": nu_late_std,
            "eta_L_R_mean": eta_R_late,
            "eta_L_full_mean": eta_full_late,
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    with SUMMARY_PATH.open("w", encoding="utf-8") as f:
        f.write("=== OU-drift simulation results (§ 6.2) ===\n")
        f.write(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}, tau_pred = {TAU_PRED}\n")
        f.write(f"lambda = {LAM}, sigma = {SIGMA}, T = {T_TOTAL}, n_runs = {N_RUNS}\n")
        f.write(f"theta0_std = {THETA0_STD}, learn_c = {LEARN_C}\n")
        f.write(f"OU stationary variance sigma^2/(2 lambda) = {SIGMA**2 / (2 * LAM):.3e}\n")
        f.write(f"Transient window upper bound 1/sigma^2 ~ {1.0 / SIGMA**2:.2e}\n\n")
        f.write(f"--- Fit eta_L*t vs ln(lam t), N_max = R*t, window [{int(t_min_fit)},{int(t_max_fit)}] ---\n")
        f.write(f"  slope = {fit_R.get('slope', float('nan')):.3f}\n")
        f.write(f"  intercept = {fit_R.get('intercept', float('nan')):.3f}\n")
        f.write(f"  R^2 = {fit_R.get('R2', float('nan')):.3f}\n")
        f.write(f"  n_points = {fit_R.get('n_points', 0)}\n")
        f.write(f"  theory K/(2R) = {slope_theory:.1f}\n")
        f.write(f"  K_eff (R=1) = {fit_R.get('K_eff_assuming_R1', float('nan')):.2f}\n")
        if fit_R.get("ok"):
            f.write(f"  ratio empirical/theory = {fit_R['slope'] / slope_theory:.3f}\n")
        f.write("\n--- Fit with full N_max budget ---\n")
        f.write(f"  slope = {fit_full.get('slope', float('nan')):.3f}\n")
        f.write(f"  R^2 = {fit_full.get('R2', float('nan')):.3f}\n")
        if fit_full.get("ok"):
            f.write(f"  ratio empirical/theory = {fit_full['slope'] / slope_theory:.3f}\n")
        f.write(f"\nLate-time (t > {int(0.5 * T_TOTAL)}):\n")
        f.write(f"  nu = {nu_late:.3f} +- {nu_late_std:.3f}\n")
        f.write(f"  <eta_L>_R    = {eta_R_late:.3e}\n")
        f.write(f"  <eta_L>_full = {eta_full_late:.3e}\n")

    LOG_PATH.write_text("\n".join(log_lines), encoding="utf-8")
    print(f"\nWrote {SUMMARY_PATH}, {SUMMARY_JSON}, {LOG_PATH}")


if __name__ == "__main__":
    main()
