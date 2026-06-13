"""Run all three experiments for paper § 6.1 and produce figs 1-3.

Usage:
    python main.py

Outputs:
    - paper/figs/fig1_eta_three_regimes.png
    - paper/figs/fig2_eta_avg_vs_tau_w.png
    - paper/figs/fig3_nu_evolution.png
    - simulations/markov_drift/results_summary.txt (numeric results for paper)
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import average_runs, simulate_run

# Resolve output paths relative to this file so the script works from any CWD.
HERE = Path(__file__).resolve().parent
PAPER_FIGS = HERE.parent.parent / "paper" / "figs"
PAPER_FIGS.mkdir(parents=True, exist_ok=True)
SUMMARY_PATH = HERE / "results_summary.txt"
SUMMARY_JSON = HERE / "results_summary.json"

# Global parameters
K = 8
TAU_PRED = 1  # one-step predictive horizon (rapid mixing of Dirichlet(alpha=0.3) chains
              # makes multi-step MI essentially zero; tau=1 is the informative scale)
ALPHA = 0.3   # Dirichlet concentration: alpha=0.3 yields tau_relax ~ 2.5 << tau_E=1000,
              # i.e. adiabatic regime lambda*tau_relax ~ 2.5e-3, consistent with § 5.2.
SEED = 20260524


def run_condition(
    label: str,
    lam: float,
    tau_w: float,
    T: int,
    n_runs: int,
    measure_every: int,
    seed_offset: int,
) -> dict:
    print(f"  [{label}] lam={lam}, tau_w={tau_w}, T={T}, runs={n_runs} ...", flush=True)
    runs = []
    for r in range(n_runs):
        rng = np.random.default_rng(SEED + seed_offset + r)
        run = simulate_run(
            k=K,
            T=T,
            lam=lam,
            tau_w=tau_w,
            tau_pred=TAU_PRED,
            rng=rng,
            measure_every=measure_every,
            alpha=ALPHA,
        )
        runs.append(run)
    avg = average_runs(runs)
    avg["raw"] = runs
    avg["label"] = label
    avg["lam"] = lam
    avg["tau_w"] = tau_w
    avg["T"] = T
    avg["n_runs"] = n_runs
    return avg


def fit_eta_t_vs_log_lam_t(times: np.ndarray, eta: np.ndarray, lam: float, t_min: float, t_max: float) -> dict:
    """Fit eta_L * t vs ln(lam * t) over [t_min, t_max]: linear slope = K/(2R).

    Returns slope, intercept, R^2, and the implied K_effective assuming R = 1.
    """
    mask = (times >= t_min) & (times <= t_max) & (eta > 0)
    if mask.sum() < 5:
        return {"ok": False, "reason": "insufficient points"}
    t = times[mask]
    y = eta[mask] * t
    x = np.log(np.maximum(lam * t, 1e-12))
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = coef
    yhat = A @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
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
    print("== Experiment A: stationary (lam=0) ==")
    cond_stat = run_condition(
        label="stationary", lam=0.0, tau_w=100, T=10_000, n_runs=80, measure_every=50, seed_offset=0
    )

    print("== Experiment B: drift, unbounded memory (lam=1e-3, tau_w=inf) ==")
    cond_drift_no_reset = run_condition(
        label="drift_no_reset", lam=1e-3, tau_w=float("inf"), T=100_000, n_runs=60,
        measure_every=200, seed_offset=1000,
    )

    print("== Experiment C: drift with reset (lam=1e-3, tau_w in {100,500,2000}) ==")
    cond_reset = {}
    for i, tau_w in enumerate([100, 500, 2000]):
        cond_reset[tau_w] = run_condition(
            label=f"drift_reset_tau_w={tau_w}", lam=1e-3, tau_w=tau_w,
            T=100_000, n_runs=60, measure_every=200, seed_offset=2000 + 1000 * i,
        )

    # Extra: tau_w sweep for fig 2
    print("== Experiment D: tau_w sweep for fig 2 ==")
    tau_w_sweep_values = [50, 100, 200, 500, 1000, 2000, 5000, 10_000]
    eta_mean_vs_tau_w = []
    eta_std_vs_tau_w = []
    for i, tau_w in enumerate(tau_w_sweep_values):
        cond = run_condition(
            label=f"sweep_tau_w={tau_w}", lam=1e-3, tau_w=tau_w, T=50_000, n_runs=40,
            measure_every=500, seed_offset=10_000 + 1000 * i,
        )
        # Take time-average of eta_L over second half (after transient)
        t = cond["times"]
        eta = cond["eta_L"]
        mask = t >= cond["T"] // 2
        # also collect per-run time averages to estimate variability across runs
        per_run = []
        for r in cond["raw"]:
            tm = r["times"] >= cond["T"] // 2
            per_run.append(float(np.mean(r["eta_L"][tm])))
        eta_mean_vs_tau_w.append(float(np.mean(per_run)))
        eta_std_vs_tau_w.append(float(np.std(per_run) / np.sqrt(len(per_run))))
        print(f"   tau_w={tau_w}: <eta_L>_late = {eta_mean_vs_tau_w[-1]:.3e} +- {eta_std_vs_tau_w[-1]:.3e}")

    # ---------- FIG 1: eta_L(t) three regimes ----------
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(cond_stat["times"], cond_stat["eta_L"], label="stationary ($\\lambda=0$, $\\tau_w=100$)",
            color="tab:blue", lw=1.5)
    ax.plot(cond_drift_no_reset["times"], cond_drift_no_reset["eta_L"],
            label="drift, no reset ($\\lambda=10^{-3}$, $\\tau_w\\to\\infty$)",
            color="tab:orange", lw=1.5)
    colors = ["tab:green", "tab:red", "tab:purple"]
    for c, tau_w in zip(colors, [100, 500, 2000]):
        d = cond_reset[tau_w]
        ax.plot(d["times"], d["eta_L"], label=f"drift + reset ($\\tau_w={tau_w}$)", color=c, lw=1.1, alpha=0.85)
    ax.set_xlabel("time $t$ (steps)")
    ax.set_ylabel(r"$\eta_L(t) = I_{\mathrm{pred}}(t,\tau) / N_{\max}(t)$")
    ax.set_title(r"Landauer efficiency $\eta_L(t)$ across three regimes (k=8, $\tau$=5)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig1_path = PAPER_FIGS / "fig1_eta_three_regimes.png"
    fig.savefig(fig1_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig1.pdf")  # vector copy for the manuscript
    plt.close(fig)
    print(f"saved {fig1_path}")

    # ---------- FIG 2: <eta_L> vs tau_w ----------
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    tau_arr = np.array(tau_w_sweep_values, dtype=float)
    eta_arr = np.array(eta_mean_vs_tau_w)
    err_arr = np.array(eta_std_vs_tau_w)
    ax.errorbar(tau_arr, eta_arr, yerr=err_arr, marker="o", capsize=3, lw=1.5, color="tab:blue")
    # mark tau_E
    tau_E = 1.0 / 1e-3
    ax.axvline(tau_E, ls="--", color="grey", alpha=0.7, label=fr"$\tau_E = 1/\lambda = {tau_E:.0f}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"memory window $\tau_w$ (steps)")
    ax.set_ylabel(r"$\langle \eta_L \rangle$ (late-time average)")
    ax.set_title(r"Time-averaged $\eta_L$ vs window $\tau_w$ at $\lambda=10^{-3}$ (k=8)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig2_path = PAPER_FIGS / "fig2_eta_avg_vs_tau_w.png"
    fig.savefig(fig2_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig2.pdf")  # vector copy for the manuscript
    plt.close(fig)
    print(f"saved {fig2_path}")

    # ---------- FIG 3: nu(t) for two regimes ----------
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(cond_drift_no_reset["times"], cond_drift_no_reset["nu"],
            label=r"drift, no reset ($\tau_w\to\infty$)", color="tab:orange", lw=1.5)
    ax.plot(cond_reset[500]["times"], cond_reset[500]["nu"],
            label=r"drift + reset ($\tau_w=500$)", color="tab:red", lw=1.1, alpha=0.85)
    ax.plot(cond_reset[2000]["times"], cond_reset[2000]["nu"],
            label=r"drift + reset ($\tau_w=2000$)", color="tab:purple", lw=1.1, alpha=0.85)
    # heuristic nu_c band
    ax.axhspan(0.5, 0.7, color="grey", alpha=0.15, label=r"heuristic $\nu_c \approx 0.6 \pm 0.1$")
    ax.set_xlabel("time $t$ (steps)")
    ax.set_ylabel(r"nostalgic fraction $\nu(t) = 1 - I_{\mathrm{pred}}/I_{\mathrm{pred}}^{\mathrm{opt}}$")
    ax.set_title(r"Evolution of nostalgic fraction $\nu(t)$ (k=8, $\lambda=10^{-3}$)")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig3_path = PAPER_FIGS / "fig3_nu_evolution.png"
    fig.savefig(fig3_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig3.pdf")  # vector copy for the manuscript
    plt.close(fig)
    print(f"saved {fig3_path}")

    # ---------- Quantitative analysis ----------
    print("\n== Quantitative analysis ==")
    # 1) Slope of eta_L * t vs ln(lam t) on drift-no-reset, transient window
    t_min = 2_000  # >> tau_E = 1000
    t_max = 50_000  # well before T
    fit_full = fit_eta_t_vs_log_lam_t(
        cond_drift_no_reset["times"], cond_drift_no_reset["eta_L"], lam=1e-3, t_min=t_min, t_max=t_max
    )
    K_theory = K * (K - 1)  # = 56
    R_used = 1.0
    slope_theory = K_theory / (2.0 * R_used)
    print(f"Fit (drift, no reset) eta_L*t vs ln(lam t) on t in [{t_min},{t_max}]:")
    print(f"  slope = {fit_full.get('slope', float('nan')):.3f}")
    print(f"  R^2 = {fit_full.get('R2', float('nan')):.3f}")
    print(f"  theory K/(2R) = {K_theory}/(2*{R_used}) = {slope_theory:.1f}")
    print(f"  ratio empirical/theory = {fit_full.get('slope', float('nan'))/slope_theory:.3f}")
    print(f"  K_eff (assuming R=1) = {fit_full.get('K_eff_assuming_R1', float('nan')):.2f} (theory: {K_theory})")

    # 2) Asymptotic nu_inf for drift-no-reset (late-time average)
    nu_drift = cond_drift_no_reset["nu"]
    t_drift = cond_drift_no_reset["times"]
    late_mask = t_drift > 60_000
    nu_inf = float(np.mean(nu_drift[late_mask])) if late_mask.any() else float("nan")
    nu_inf_std = float(np.std(nu_drift[late_mask])) if late_mask.any() else float("nan")
    print(f"\nNu_inf (drift, no reset, late-time mean t>60k): {nu_inf:.3f} +- {nu_inf_std:.3f}")

    # 3) Find tau_w* maximizing <eta_L>
    best_idx = int(np.argmax(eta_arr))
    print(f"\nOptimal tau_w in sweep: {tau_w_sweep_values[best_idx]} (eta = {eta_arr[best_idx]:.3e})")

    # 4) Stationary eta_L plateau
    t_stat = cond_stat["times"]
    plateau_mask = t_stat > 5_000
    eta_plateau = float(np.mean(cond_stat["eta_L"][plateau_mask]))
    print(f"\nStationary eta_L late-time plateau: {eta_plateau:.3e}")

    # 5) Critical nu where eta_L crosses below half its peak (drift no reset)
    eta_no_reset = cond_drift_no_reset["eta_L"]
    nu_no_reset = cond_drift_no_reset["nu"]
    if eta_no_reset.max() > 0:
        i_peak = int(np.argmax(eta_no_reset))
        nu_at_peak = float(nu_no_reset[i_peak])
        half_target = eta_no_reset[i_peak] / 2.0
        post_peak = eta_no_reset[i_peak:]
        below = np.where(post_peak < half_target)[0]
        if below.size:
            j = int(below[0]) + i_peak
            nu_at_halflife = float(nu_no_reset[j])
        else:
            nu_at_halflife = float("nan")
    else:
        nu_at_peak = nu_at_halflife = float("nan")
    print(f"\nNu at eta_L peak (drift no reset): {nu_at_peak:.3f}")
    print(f"Nu when eta_L halves from peak: {nu_at_halflife:.3f}  <-- empirical nu_c proxy")

    # ---------- Write summary text and json ----------
    summary = {
        "K": K,
        "K_theory_k_times_km1": K_theory,
        "tau_pred": TAU_PRED,
        "stationary_eta_plateau": eta_plateau,
        "drift_no_reset": {
            "fit_eta_t_vs_log_lam_t": fit_full,
            "slope_theory_K_over_2R": slope_theory,
            "nu_inf": nu_inf,
            "nu_inf_std": nu_inf_std,
            "nu_at_peak": nu_at_peak,
            "nu_at_eta_halflife": nu_at_halflife,
        },
        "tau_w_sweep": {
            "tau_w_values": tau_w_sweep_values,
            "eta_mean_late": eta_mean_vs_tau_w,
            "eta_sem_late": eta_std_vs_tau_w,
            "optimal_tau_w": tau_w_sweep_values[best_idx],
            "optimal_eta": float(eta_arr[best_idx]),
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    with SUMMARY_PATH.open("w", encoding="utf-8") as f:
        f.write("=== Markov-drift simulation results (§ 6.1) ===\n")
        f.write(f"k = {K}, K = k(k-1) = {K_theory}, tau_pred = {TAU_PRED}\n")
        f.write(f"Stationary eta_L late-time plateau: {eta_plateau:.3e}\n\n")
        f.write("Drift, no reset (lambda=1e-3, tau_w=inf, T=1e5, 60 runs):\n")
        f.write(f"  fit eta_L*t vs ln(lam t) on t in [{t_min},{t_max}]:\n")
        f.write(f"    slope = {fit_full.get('slope', float('nan')):.3f}\n")
        f.write(f"    R^2 = {fit_full.get('R2', float('nan')):.3f}\n")
        f.write(f"    theory K/(2R) = {slope_theory:.1f}\n")
        f.write(f"    K_eff (R=1) = {fit_full.get('K_eff_assuming_R1', float('nan')):.2f}\n")
        f.write(f"  nu_inf (late-time t>60k): {nu_inf:.3f} +- {nu_inf_std:.3f}\n")
        f.write(f"  nu at eta_L peak: {nu_at_peak:.3f}\n")
        f.write(f"  nu at eta_L halflife (empirical nu_c proxy): {nu_at_halflife:.3f}\n\n")
        f.write("Tau_w sweep (lambda=1e-3, T=5e4, 40 runs):\n")
        for tw, em, es in zip(tau_w_sweep_values, eta_mean_vs_tau_w, eta_std_vs_tau_w):
            f.write(f"  tau_w={tw:>6}: <eta_L>_late = {em:.3e} +- {es:.3e}\n")
        f.write(f"  optimal tau_w in grid: {tau_w_sweep_values[best_idx]}\n")
    print(f"\nWrote {SUMMARY_PATH} and {SUMMARY_JSON}")


if __name__ == "__main__":
    main()
