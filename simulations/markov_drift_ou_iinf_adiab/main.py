"""iter-006 Step 2: ADIABATIC SCAN of OU + I_infinity (BNT redundancy).

Background (from Step 1, ../markov_drift_ou_iinf/results_summary.txt):
    The joint fit  L_excess(t) = A*t + B*ln(lambda*t) + C  on Step 1 data
    (sigma=0.01, sigma^2/lambda=0.1) gave A=1.85e-2, B=9.5, R^2=0.9998.
    The functional form is excellent but B is only ~34% of the BNT-predicted
    K/2 = 28 for K=k(k-1)=56 with k=8.

Hypothesis: in OU drift the cumulative KL has two contributions,
   (a) a linear drift-loss A*t from structural unlearnability of the moving
       parameter (proportional to sigma^2/lambda in some regime),
   (b) the BNT log-redundancy term (K/2) ln(lambda*t) from posterior
       concentration on the transient time-scale 1/sigma^2.

If (a) "drowns" (b) at sigma=0.01, then driving sigma down (adiabatic limit
sigma^2/lambda -> 0) should let B saturate at K/2.

This script scans  sigma in {0.01, 0.003, 0.001} (sigma^2/lambda decreasing by
two orders of magnitude) and measures B at each setting.

NOTE on the result of this scan (see results_summary.txt): the joint fit on a
fixed comparison window [1e3, 1e4] does show B drifting upward as sigma
decreases, but it does NOT settle at K/2 = 28; instead it keeps growing and
crosses through K/2 on its way to larger values when the fit window is widened
toward 0.5/sigma^2. This is interpreted in the summary.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import simulate_batch_iinf

HERE = Path(__file__).resolve().parent
SUMMARY_PATH = HERE / "results_summary.txt"
SUMMARY_JSON = HERE / "results_summary.json"
LOG_PATH = HERE / "run.log"

# ---------------- Fixed parameters ----------------
K_STATES = 8
K_PARAMS = K_STATES * (K_STATES - 1)
LAM = 1e-3
N_RUNS = 25
THETA0_STD = 0.1
LEARN_C = 3.0
R = 1.0
SEED = 20260524

# ---------------- Scan settings ----------------
# (sigma, T_total, measure_every, label)
SCAN = [
    (0.01,  50_000,  100, "control_iter006_step1_reproduction"),
    (0.003, 50_000,  100, "intermediate"),
    (0.001, 200_000, 200, "adiabatic"),
]

# Comparison fit window (matches Step 1): always [1e3, 1e4]
T_MIN_COMP = 1_000.0
T_MAX_COMP = 10_000.0


# --------------------------------------------------------------------------- #
# Fit helpers                                                                 #
# --------------------------------------------------------------------------- #

def joint_fit(times: np.ndarray, y: np.ndarray, lam: float,
              t_min: float, t_max: float) -> dict:
    """Fit y(t) = A*t + B*ln(lam*t) + C on [t_min, t_max]."""
    mask = (times >= t_min) & (times <= t_max)
    t = times[mask]
    if t.size < 5:
        return {"ok": False, "n_points": int(t.size),
                "t_min": float(t_min), "t_max": float(t_max)}
    yy = y[mask]
    ln = np.log(np.maximum(lam * t, 1e-12))
    A_mat = np.vstack([t, ln, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A_mat, yy, rcond=None)
    yhat = A_mat @ coef
    ss_res = float(np.sum((yy - yhat) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "ok": True, "n_points": int(t.size),
        "A_linear": float(coef[0]), "B_log": float(coef[1]), "C": float(coef[2]),
        "R2": float(r2), "t_min": float(t_min), "t_max": float(t_max),
    }


def log_only_fit(times: np.ndarray, y: np.ndarray, lam: float,
                 t_min: float, t_max: float) -> dict:
    """Fit y(t) = B*ln(lam*t) + C on [t_min, t_max]."""
    mask = (times >= t_min) & (times <= t_max)
    t = times[mask]
    if t.size < 5:
        return {"ok": False, "n_points": int(t.size)}
    yy = y[mask]
    ln = np.log(np.maximum(lam * t, 1e-12))
    A_mat = np.vstack([ln, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A_mat, yy, rcond=None)
    yhat = A_mat @ coef
    ss_res = float(np.sum((yy - yhat) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "ok": True, "n_points": int(t.size),
        "B_log": float(coef[0]), "C": float(coef[1]),
        "R2": float(r2), "t_min": float(t_min), "t_max": float(t_max),
    }


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #

def main() -> None:
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_lines.append(msg)

    log("== iter-006 Step 2: OU + I_infinity ADIABATIC SCAN ==")
    log(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}, K/2 = {K_PARAMS/2}")
    log(f"lambda = {LAM}, theta0_std = {THETA0_STD}, learn_c = {LEARN_C}")
    log(f"N_RUNS = {N_RUNS}, seed = {SEED}")
    log(f"Scan settings: {[(s, T, me, lbl) for s, T, me, lbl in SCAN]}")
    log(f"Comparison fit window (matches Step 1): [{T_MIN_COMP:.0f}, {T_MAX_COMP:.0f}]")

    t_start = time.time()
    results = []
    for sigma, T_total, measure_every, label in SCAN:
        sigma2_over_lambda = sigma ** 2 / LAM
        transient_bound = 1.0 / sigma ** 2

        log(f"\n--- Running sigma = {sigma:g}  ({label}) ---")
        log(f"   sigma^2/lambda = {sigma2_over_lambda:.3e}, "
            f"transient T_max = 1/sigma^2 = {transient_bound:.2e}")
        log(f"   T_total = {T_total}, measure_every = {measure_every}, N_runs = {N_RUNS}")

        rng = np.random.default_rng(SEED + int(round(1000 * sigma)))
        t0 = time.time()
        out = simulate_batch_iinf(
            k=K_STATES, T=T_total, lam=LAM, sigma=sigma, n_runs=N_RUNS, rng=rng,
            measure_every=measure_every, theta0_std=THETA0_STD,
            learn_c=LEARN_C, R=R,
        )
        dt = time.time() - t0
        log(f"   simulated in {dt:.1f}s ({T_total*N_RUNS/dt/1e6:.2f} Msteps/s)")

        times = out["times"]
        L_excess = out["L_excess_cum"]
        L_excess_std = out["L_excess_cum_std"]

        # Fit window 1: comparison window [1e3, 1e4] (matches Step 1)
        # This is fully *inside* the transient bound for all sigma in our scan.
        fit_comp = joint_fit(times, L_excess, LAM, T_MIN_COMP, T_MAX_COMP)
        # Fit window 2: principled "as wide as allowed" window
        # [1e3, min(T_total, 0.5/sigma^2)]
        t_max_wide = min(float(T_total), 0.5 * transient_bound)
        fit_wide = joint_fit(times, L_excess, LAM, T_MIN_COMP, t_max_wide)
        # Auxiliary: log-only fit on the comparison window
        fit_log_comp = log_only_fit(times, L_excess, LAM, T_MIN_COMP, T_MAX_COMP)

        log(f"   joint fit (comparison [{T_MIN_COMP:.0f},{T_MAX_COMP:.0f}]):")
        log(f"      A = {fit_comp['A_linear']:.4e},  B = {fit_comp['B_log']:.3f},  "
            f"C = {fit_comp['C']:.3f},  R^2 = {fit_comp['R2']:.4f},  "
            f"B/(K/2) = {fit_comp['B_log']/(K_PARAMS/2):.3f}")
        log(f"   joint fit (wide [{T_MIN_COMP:.0f},{t_max_wide:.0f}]):")
        log(f"      A = {fit_wide['A_linear']:.4e},  B = {fit_wide['B_log']:.3f},  "
            f"C = {fit_wide['C']:.3f},  R^2 = {fit_wide['R2']:.4f},  "
            f"B/(K/2) = {fit_wide['B_log']/(K_PARAMS/2):.3f}")
        log(f"   log-only fit (no A*t term, comparison window):  B = {fit_log_comp['B_log']:.3f}, "
            f"R^2 = {fit_log_comp['R2']:.4f}")

        # Diagnostic: drift-loss dominance ratio A*t / L_excess at several t
        diag_ratio = {}
        for t_check in [1e3, 1e4, 5e4, 1e5, 2e5]:
            if t_check > times[-1]: continue
            i = int(np.argmin(np.abs(times - t_check)))
            A_t = fit_wide["A_linear"] * times[i]
            tot = L_excess[i] if L_excess[i] != 0 else 1.0
            diag_ratio[int(times[i])] = float(A_t / tot)

        results.append({
            "sigma": sigma,
            "label": label,
            "sigma2_over_lambda": sigma2_over_lambda,
            "transient_bound": transient_bound,
            "T_total": T_total,
            "measure_every": measure_every,
            "runtime_s": dt,
            "fit_comparison_window": fit_comp,
            "fit_wide_window": fit_wide,
            "fit_log_only_comp_window": fit_log_comp,
            "drift_loss_ratio_A_times_t_over_L": diag_ratio,
            "times": times.tolist(),
            "L_excess_cum": L_excess.tolist(),
            "L_excess_cum_std": L_excess_std.tolist(),
            "I_obs_1step": out["I_obs_1step"].tolist(),
            "I_opt_1step": out["I_opt_1step"].tolist(),
            "nu": out["nu"].tolist(),
            "theta_err": out["theta_err"].tolist(),
        })

    total_runtime = time.time() - t_start
    log(f"\n== Total scan runtime: {total_runtime:.1f}s ==")

    # ---------------- Plot 1: B vs sigma^2/lambda ----------------
    sigma2_lambdas = np.array([r["sigma2_over_lambda"] for r in results])
    B_comp = np.array([r["fit_comparison_window"]["B_log"] for r in results])
    B_wide = np.array([r["fit_wide_window"]["B_log"] for r in results])
    B_logonly = np.array([r["fit_log_only_comp_window"]["B_log"] for r in results])
    sigmas = np.array([r["sigma"] for r in results])

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.semilogx(sigma2_lambdas, B_comp, "o-", color="tab:blue", ms=8, lw=1.5,
                label=rf"joint fit $A\cdot t + B\ln(\lambda t)$, comparison window "
                      rf"$[10^3, 10^4]$")
    ax.semilogx(sigma2_lambdas, B_wide, "s--", color="tab:orange", ms=7, lw=1.2,
                label=rf"joint fit, wide window $[10^3, 0.5/\sigma^2]$")
    ax.semilogx(sigma2_lambdas, B_logonly, "v:", color="tab:green", ms=6, lw=1.0,
                label=r"log-only fit $B\ln(\lambda t)$, comparison window")
    ax.axhline(K_PARAMS / 2, color="tab:red", ls="--", lw=1.5,
               label=rf"BNT theory $K/2 = {K_PARAMS/2:.0f}$")
    ax.set_xlabel(r"$\sigma^2/\lambda$  (adiabatic parameter, $\to 0$ = adiabatic)")
    ax.set_ylabel(r"fitted log coefficient $B$")
    ax.set_title(r"Adiabatic scan: BNT log-coefficient $B$ vs $\sigma^2/\lambda$ "
                 rf"($N_{{\rm runs}}={N_RUNS}$)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    # annotate sigma values
    for sig, x, y in zip(sigmas, sigma2_lambdas, B_comp):
        ax.annotate(rf"$\sigma={sig:g}$", (x, y),
                    textcoords="offset points", xytext=(5, 6), fontsize=8)
    fig.tight_layout()
    fig1_path = HERE / "fig_adiab_B_vs_sigma2_lambda.png"
    fig.savefig(fig1_path, dpi=150)
    fig.savefig(HERE.parent.parent / "paper" / "figs" / "Fig6.pdf")  # vector copy for the supplementary
    plt.close(fig)
    log(f"saved {fig1_path}")

    # ---------------- Plot 2: L_excess curves overlaid + drift-loss ratio ----------------
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0))

    ax = axes[0]
    colours = ["tab:blue", "tab:orange", "tab:green"]
    for r, colour in zip(results, colours):
        times = np.array(r["times"])
        L = np.array(r["L_excess_cum"])
        ax.plot(np.log(LAM * times), L, "-", color=colour, lw=1.5,
                label=rf"$\sigma={r['sigma']:g}$ ($\sigma^2/\lambda={r['sigma2_over_lambda']:.1e}$)")
        # overlay the joint fit on the wide window
        f = r["fit_wide_window"]
        if f.get("ok"):
            mask = (times >= f["t_min"]) & (times <= f["t_max"])
            tt = times[mask]
            y_th = f["A_linear"] * tt + f["B_log"] * np.log(LAM * tt) + f["C"]
            ax.plot(np.log(LAM * tt), y_th, "--", color=colour, lw=0.8, alpha=0.7)
    ax.set_xlabel(r"$\ln(\lambda t)$")
    ax.set_ylabel(r"$L_{\rm excess}(t)$  [nats]")
    ax.set_title(r"$L_{\rm excess}(t)$ for adiabatic scan (dashed = joint fit on wide window)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for r, colour in zip(results, colours):
        times = np.array(r["times"])
        L = np.array(r["L_excess_cum"])
        A = r["fit_wide_window"]["A_linear"]
        ratio = (A * times) / np.maximum(L, 1e-9)
        ax.plot(times, ratio, "-", color=colour, lw=1.4,
                label=rf"$\sigma={r['sigma']:g}$, $A={A:.2e}$")
    ax.axhline(1.0, color="grey", ls=":", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$A\cdot t / L_{\rm excess}(t)$  (drift-loss share of cumulative KL)")
    ax.set_title("Drift-loss vs BNT-log decomposition: ratio of $A\\cdot t$ to total $L_{\\rm excess}$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig2_path = HERE / "fig_adiab_Lexcess_scan.png"
    fig.savefig(fig2_path, dpi=150)
    fig.savefig(HERE.parent.parent / "paper" / "figs" / "Fig7.pdf")  # vector copy for the supplementary
    plt.close(fig)
    log(f"saved {fig2_path}")

    # ---------------- Summary outputs ----------------
    summary = {
        "experiment": "iter-006 step 2: OU + I_infinity adiabatic scan",
        "k": K_STATES,
        "K_params": K_PARAMS,
        "K_over_2": K_PARAMS / 2,
        "lambda": LAM,
        "n_runs": N_RUNS,
        "theta0_std": THETA0_STD,
        "learn_c": LEARN_C,
        "seed": SEED,
        "total_runtime_s": total_runtime,
        "comparison_fit_window": [T_MIN_COMP, T_MAX_COMP],
        "scan": [
            {**{kk: vv for kk, vv in r.items()
                if kk not in ("times", "L_excess_cum", "L_excess_cum_std",
                              "I_obs_1step", "I_opt_1step", "nu", "theta_err")}}
            for r in results
        ],
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    with SUMMARY_PATH.open("w", encoding="utf-8") as f:
        f.write("=== iter-006 Step 2: OU + I_infinity ADIABATIC SCAN ===\n")
        f.write(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}, K/2 = {K_PARAMS/2}\n")
        f.write(f"lambda = {LAM}, N_runs = {N_RUNS}, learn_c = {LEARN_C}, seed = {SEED}\n")
        f.write(f"Total runtime: {total_runtime:.1f}s\n\n")

        f.write("--- Adiabatic scan summary ---\n")
        f.write(f"{'sigma':>7s}  {'sig^2/lam':>10s}  {'1/sig^2':>11s}  "
                f"{'T_tot':>7s}  {'B_comp':>7s}  {'B_wide':>7s}  "
                f"{'B/(K/2)_c':>9s}  {'B/(K/2)_w':>9s}  "
                f"{'R2_comp':>7s}  {'R2_wide':>7s}\n")
        for r in results:
            fc = r["fit_comparison_window"]
            fw = r["fit_wide_window"]
            f.write(f"{r['sigma']:>7g}  {r['sigma2_over_lambda']:>10.3e}  "
                    f"{r['transient_bound']:>11.2e}  {r['T_total']:>7d}  "
                    f"{fc['B_log']:>7.3f}  {fw['B_log']:>7.3f}  "
                    f"{fc['B_log']/(K_PARAMS/2):>9.3f}  {fw['B_log']/(K_PARAMS/2):>9.3f}  "
                    f"{fc['R2']:>7.4f}  {fw['R2']:>7.4f}\n")

        f.write("\nFull fit details:\n")
        for r in results:
            f.write(f"\n  sigma = {r['sigma']:g} ({r['label']})\n")
            f.write(f"    sigma^2/lambda = {r['sigma2_over_lambda']:.3e}\n")
            f.write(f"    transient bound 1/sigma^2 = {r['transient_bound']:.2e}\n")
            f.write(f"    T_total = {r['T_total']}, measure_every = {r['measure_every']}\n")
            f.write(f"    runtime = {r['runtime_s']:.1f}s\n")
            fc = r["fit_comparison_window"]
            f.write(f"    joint fit [comparison {fc['t_min']:.0f},{fc['t_max']:.0f}]:\n")
            f.write(f"       A_linear = {fc['A_linear']:.5e}\n")
            f.write(f"       B_log    = {fc['B_log']:.3f}\n")
            f.write(f"       C        = {fc['C']:.3f}\n")
            f.write(f"       R^2      = {fc['R2']:.4f}, n_pts = {fc['n_points']}\n")
            f.write(f"       B/(K/2)  = {fc['B_log']/(K_PARAMS/2):.3f}\n")
            fw = r["fit_wide_window"]
            f.write(f"    joint fit [wide       {fw['t_min']:.0f},{fw['t_max']:.0f}]:\n")
            f.write(f"       A_linear = {fw['A_linear']:.5e}\n")
            f.write(f"       B_log    = {fw['B_log']:.3f}\n")
            f.write(f"       C        = {fw['C']:.3f}\n")
            f.write(f"       R^2      = {fw['R2']:.4f}, n_pts = {fw['n_points']}\n")
            f.write(f"       B/(K/2)  = {fw['B_log']/(K_PARAMS/2):.3f}\n")
            fl = r["fit_log_only_comp_window"]
            f.write(f"    log-only fit [comparison]:\n")
            f.write(f"       B_log    = {fl['B_log']:.3f}, R^2 = {fl['R2']:.4f}\n")
            f.write(f"    drift-loss ratio A*t / L_excess(t):\n")
            for t_c, rat in r["drift_loss_ratio_A_times_t_over_L"].items():
                f.write(f"       t={t_c:>7d}:  ratio = {rat:.3f}\n")

        f.write("\n--- INTERPRETATION ---\n")
        f.write("See README.md for the qualitative conclusion driving the article-level\n")
        f.write("decision (does conjecture S8.1 / (S8.2) Supplementary survive as adiabatic-only\n")
        f.write("or does it need fallback).\n")

    LOG_PATH.write_text("\n".join(log_lines), encoding="utf-8")
    print(f"\nWrote {SUMMARY_PATH}, {SUMMARY_JSON}, {LOG_PATH}")


if __name__ == "__main__":
    main()
