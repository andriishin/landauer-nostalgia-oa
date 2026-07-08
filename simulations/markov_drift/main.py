"""Run all three experiments for paper § 6.1 and produce figs 1-3.

Re-equipped apparatus: efficiency eta = I_pred / I_mem (was I_pred / N_max).
The old "slope K/(2R)" fit was an artefact of N_max = R*t and has been removed.

Usage:
    python main.py

Outputs:
    - paper/figs/fig1_Ipred_three_regimes.png  (+ Fig1.pdf)
    - paper/figs/fig2_Ipred_avg_vs_tau_w.png   (+ Fig2.pdf)
    - paper/figs/fig3_nu_evolution.png         (+ Fig3.pdf)
    - simulations/markov_drift/results_summary.txt / .json
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
K_PARAMS = K * (K - 1)  # = 56 free off-diagonal transition parameters
TAU_PRED = 1  # one-step predictive horizon (rapid mixing of Dirichlet(alpha=0.3) chains
              # makes multi-step MI essentially zero; tau=1 is the informative scale)
ALPHA = 0.3   # Dirichlet concentration: alpha=0.3 yields tau_relax ~ 2.5 << tau_E=1000,
              # i.e. adiabatic regime lambda*tau_relax ~ 2.5e-3, consistent with § 5.2.
LAM = 1e-3
TAU_E = 1.0 / LAM  # environmental coherence time = 1/lambda = 1000 steps
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


def late_time_mean(cond: dict, key: str, frac: float = 0.5) -> tuple[float, float]:
    """Per-run late-time average of `key` over the last `frac` of the run, then
    mean +- SEM across runs. Returns (mean, sem)."""
    T = cond["T"]
    per_run = []
    for r in cond["raw"]:
        m = r["times"] >= T * (1.0 - frac)
        if m.any():
            per_run.append(float(np.mean(r[key][m])))
    if not per_run:
        return float("nan"), float("nan")
    arr = np.array(per_run)
    return float(arr.mean()), float(arr.std() / np.sqrt(len(arr)))


def main() -> None:
    print("== Experiment A: stationary (lam=0) ==")
    cond_stat = run_condition(
        label="stationary", lam=0.0, tau_w=float("inf"), T=10_000, n_runs=80,
        measure_every=50, seed_offset=0,
    )

    print("== Experiment B: drift, unbounded memory (lam=1e-3, tau_w=inf) ==")
    cond_drift_no_reset = run_condition(
        label="drift_no_reset", lam=LAM, tau_w=float("inf"), T=100_000, n_runs=60,
        measure_every=200, seed_offset=1000,
    )

    print("== Experiment C: drift with reset (lam=1e-3, tau_w in {100,500,2000}) ==")
    cond_reset = {}
    for i, tau_w in enumerate([100, 500, 2000]):
        cond_reset[tau_w] = run_condition(
            label=f"drift_reset_tau_w={tau_w}", lam=LAM, tau_w=tau_w,
            T=100_000, n_runs=60, measure_every=200, seed_offset=2000 + 1000 * i,
        )

    # Extra: tau_w sweep for fig 2 (optimal forgetting)
    print("== Experiment D: tau_w sweep for fig 2 ==")
    tau_w_sweep_values = [50, 100, 200, 500, 1000, 2000, 5000, 10_000]
    Ipred_mean_vs_tau_w, Ipred_sem_vs_tau_w = [], []
    eta_mean_vs_tau_w, eta_sem_vs_tau_w = [], []
    for i, tau_w in enumerate(tau_w_sweep_values):
        cond = run_condition(
            label=f"sweep_tau_w={tau_w}", lam=LAM, tau_w=tau_w, T=50_000, n_runs=40,
            measure_every=500, seed_offset=10_000 + 1000 * i,
        )
        ip_m, ip_s = late_time_mean(cond, "I_pred")
        et_m, et_s = late_time_mean(cond, "eta")
        Ipred_mean_vs_tau_w.append(ip_m)
        Ipred_sem_vs_tau_w.append(ip_s)
        eta_mean_vs_tau_w.append(et_m)
        eta_sem_vs_tau_w.append(et_s)
        print(f"   tau_w={tau_w:>6}: <I_pred>_late = {ip_m:.3e} +- {ip_s:.1e} | "
              f"<eta>_late = {et_m:.3e} +- {et_s:.1e}")

    # ---------- FIG 1: I_pred(t) three regimes ----------
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(cond_stat["times"], cond_stat["I_pred"],
            label=r"stationary ($\lambda=0$, $\tau_w\to\infty$)", color="tab:blue", lw=1.5)
    ax.plot(cond_drift_no_reset["times"], cond_drift_no_reset["I_pred"],
            label=r"drift, no reset ($\lambda=10^{-3}$, $\tau_w\to\infty$)",
            color="tab:orange", lw=1.5)
    colors = ["tab:green", "tab:red", "tab:purple"]
    for c, tau_w in zip(colors, [100, 500, 2000]):
        d = cond_reset[tau_w]
        ax.plot(d["times"], d["I_pred"], label=fr"drift + reset ($\tau_w={tau_w}$)",
                color=c, lw=1.1, alpha=0.85)
    ax.set_xlabel("time $t$ (steps)")
    ax.set_ylabel(r"predictive information $I_{\mathrm{pred}}(t,\tau)$ (nats)")
    ax.set_title(r"Collapse vs retention of $I_{\mathrm{pred}}(t)$ "
                 r"(k=8, $\tau$=1, $\lambda=10^{-3}$)")
    ax.set_xscale("log")
    ax.axvline(TAU_E, ls=":", color="grey", alpha=0.6)
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig1_path = PAPER_FIGS / "fig1_Ipred_three_regimes.png"
    fig.savefig(fig1_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig1.pdf")
    plt.close(fig)
    print(f"saved {fig1_path}")

    # ---------- FIG 2: optimal forgetting -- <I_pred> (and <eta>) vs tau_w ----------
    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    tau_arr = np.array(tau_w_sweep_values, dtype=float)
    ip_arr = np.array(Ipred_mean_vs_tau_w)
    ip_err = np.array(Ipred_sem_vs_tau_w)
    et_arr = np.array(eta_mean_vs_tau_w)
    et_err = np.array(eta_sem_vs_tau_w)

    ax.errorbar(tau_arr, ip_arr, yerr=ip_err, marker="o", capsize=3, lw=1.5,
                color="tab:blue", label=r"$\langle I_{\mathrm{pred}}\rangle$ (nats)")
    ax.set_xscale("log")
    ax.set_xlabel(r"memory window $\tau_w$ (steps)")
    ax.set_ylabel(r"$\langle I_{\mathrm{pred}}\rangle$ (late-time)", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax.twinx()
    ax2.errorbar(tau_arr, et_arr, yerr=et_err, marker="s", capsize=3, lw=1.2,
                 color="tab:green", alpha=0.8, label=r"$\langle\eta\rangle$")
    ax2.set_ylabel(r"$\langle\eta\rangle = \langle I_{\mathrm{pred}}/I_{\mathrm{mem}}\rangle$",
                   color="tab:green")
    ax2.tick_params(axis="y", labelcolor="tab:green")

    # mark optimum (by I_pred) and tau_E
    best_idx = int(np.argmax(ip_arr))
    ax.axvline(tau_arr[best_idx], ls="-", color="tab:blue", alpha=0.35,
               label=fr"$\tau_w^* = {tau_w_sweep_values[best_idx]}$")
    ax.axvline(TAU_E, ls="--", color="grey", alpha=0.7,
               label=fr"$\tau_E = 1/\lambda = {TAU_E:.0f}$")
    ax.set_title(r"Optimal forgetting: $\langle I_{\mathrm{pred}}\rangle$ vs $\tau_w$ "
                 r"at $\lambda=10^{-3}$ (k=8)")
    ax.grid(True, which="both", alpha=0.3)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=8)
    fig.tight_layout()
    fig2_path = PAPER_FIGS / "fig2_Ipred_avg_vs_tau_w.png"
    fig.savefig(fig2_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig2.pdf")
    plt.close(fig)
    print(f"saved {fig2_path}")

    # ---------- FIG 3: nu_Still(t) and nu_op(t) ----------
    # Memory-ballast floor c = M_refresh * tau_E / |M| (Lemma 2).
    # Refresh rate M_refresh: a sliding window of length tau_w refreshes its whole
    # content over tau_w steps, so it ingests ~1 transition/step; over the coherence
    # time tau_E it refreshes min(tau_w, tau_E) of its |M| ~ tau_w slots. We report the
    # information-side estimate c ~ 1 - <eta> directly from data (see summary).
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(cond_drift_no_reset["times"], cond_drift_no_reset["nu_still"],
            label=r"$\nu^{\mathrm{Still}}$, drift, no reset ($\tau_w\to\infty$)",
            color="tab:orange", lw=1.6)
    ax.plot(cond_reset[500]["times"], cond_reset[500]["nu_still"],
            label=r"$\nu^{\mathrm{Still}}$, drift + reset ($\tau_w=500$)",
            color="tab:red", lw=1.2, alpha=0.9)
    ax.plot(cond_reset[2000]["times"], cond_reset[2000]["nu_still"],
            label=r"$\nu^{\mathrm{Still}}$, drift + reset ($\tau_w=2000$)",
            color="tab:purple", lw=1.2, alpha=0.9)
    # secondary diagnostic nu_op (dashed)
    ax.plot(cond_drift_no_reset["times"], cond_drift_no_reset["nu_op"],
            label=r"$\nu^{\mathrm{op}}$, drift, no reset", color="tab:orange",
            lw=1.0, ls="--", alpha=0.7)
    ax.plot(cond_reset[500]["times"], cond_reset[500]["nu_op"],
            label=r"$\nu^{\mathrm{op}}$, drift + reset ($\tau_w=500$)",
            color="tab:red", lw=1.0, ls="--", alpha=0.7)
    ax.axhline(1.0, ls=":", color="grey", alpha=0.6)
    ax.set_xlabel("time $t$ (steps)")
    ax.set_ylabel(r"nostalgic fractions $\nu^{\mathrm{Still}}=1-\eta$, "
                  r"$\nu^{\mathrm{op}}=1-I_{\mathrm{pred}}/I_{\mathrm{opt}}$")
    ax.set_title(r"Nostalgia: $\nu^{\mathrm{Still}}\to 1$ under unbounded memory, "
                 r"bounded under reset (k=8, $\lambda=10^{-3}$)")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="best", fontsize=7.5)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig3_path = PAPER_FIGS / "fig3_nu_evolution.png"
    fig.savefig(fig3_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig3.pdf")
    plt.close(fig)
    print(f"saved {fig3_path}")

    # ---------- Quantitative analysis ----------
    print("\n== Quantitative analysis (late-time, last 40%) ==")

    def regime_block(cond: dict, frac: float = 0.4) -> dict:
        return {
            "I_pred": late_time_mean(cond, "I_pred", frac),
            "I_mem": late_time_mean(cond, "I_mem", frac),
            "eta": late_time_mean(cond, "eta", frac),
            "nu_still": late_time_mean(cond, "nu_still", frac),
            "nu_op": late_time_mean(cond, "nu_op", frac),
            "n_obs": late_time_mean(cond, "n_obs", frac),
        }

    blocks = {
        "stationary": regime_block(cond_stat),
        "drift_no_reset": regime_block(cond_drift_no_reset),
        "drift_reset_100": regime_block(cond_reset[100]),
        "drift_reset_500": regime_block(cond_reset[500]),
        "drift_reset_2000": regime_block(cond_reset[2000]),
    }
    for name, b in blocks.items():
        print(f"  [{name}]")
        for q in ("I_pred", "I_mem", "eta", "nu_still", "nu_op"):
            m, s = b[q]
            print(f"      <{q}> = {m:.4f} +- {s:.4f}")

    # I_pred collapse check: drift-no-reset late I_pred vs stationary I_pred
    ip_stat = blocks["stationary"]["I_pred"][0]
    ip_drift = blocks["drift_no_reset"]["I_pred"][0]
    collapse_ratio = ip_drift / ip_stat if ip_stat > 0 else float("nan")
    print(f"\n  I_pred collapse: drift-no-reset/stationary = {collapse_ratio:.4f} "
          f"({ip_drift:.4f} / {ip_stat:.4f})")

    # DPI / eta-in-[0,1] check across every sampled point of every condition
    all_conds = [cond_stat, cond_drift_no_reset] + list(cond_reset.values())
    eta_min, eta_max, dpi_violations = np.inf, -np.inf, 0
    total_points = 0
    for cond in all_conds:
        for r in cond["raw"]:
            e = r["eta"]
            eta_min = min(eta_min, float(e.min()))
            eta_max = max(eta_max, float(e.max()))
            dpi_violations += int(np.sum(r["I_pred"] > r["I_mem"] + 1e-9))
            total_points += e.size
    print(f"\n  eta range over ALL {total_points} sampled points: "
          f"[{eta_min:.4f}, {eta_max:.4f}]")
    print(f"  DPI violations (I_pred > I_mem): {dpi_violations}")

    # Optimal tau_w (by late <I_pred>) and by late <eta>
    best_ip_idx = int(np.argmax(np.array(Ipred_mean_vs_tau_w)))
    best_eta_idx = int(np.argmax(np.array(eta_mean_vs_tau_w)))
    tau_w_star_Ipred = tau_w_sweep_values[best_ip_idx]
    tau_w_star_eta = tau_w_sweep_values[best_eta_idx]
    print(f"\n  Optimal tau_w by <I_pred>: {tau_w_star_Ipred} "
          f"(tau_E = {TAU_E:.0f})")
    print(f"  Optimal tau_w by <eta>:    {tau_w_star_eta}")

    # Memory-ballast floor c (Lemma 2) and check nu >= 1 - c.
    # Information-side estimate: at the optimum the retained model achieves
    # <eta>* = <I_pred/I_mem>*, so the *minimum attainable* ballast is
    # c_floor = 1 - max_tau_w <eta>; any unbounded-memory run must obey
    # nu_Still >= this floor (it cannot do better than the best window).
    eta_star = float(np.max(np.array(eta_mean_vs_tau_w)))
    c_floor = 1.0 - eta_star
    nu_still_drift = blocks["drift_no_reset"]["nu_still"][0]
    floor_ok = nu_still_drift >= c_floor - 1e-6
    print(f"\n  Lemma 2 floor: c = 1 - max<eta> = 1 - {eta_star:.4f} = {c_floor:.4f}")
    print(f"  drift-no-reset <nu_Still> = {nu_still_drift:.4f} "
          f">= 1 - c = {c_floor:.4f} ? {floor_ok}")

    # ---------- Write summary text and json ----------
    def pair(b, q):
        return {"mean": b[q][0], "sem": b[q][1]}

    summary = {
        "apparatus": "eta = I_pred / I_mem ;  I_mem = (K/2) ln(max(N_obs,e)) ;  "
                     "nu_Still = 1 - eta ;  nu_op = 1 - I_pred/I_opt",
        "k": K,
        "K_params_k_times_km1": K_PARAMS,
        "tau_pred": TAU_PRED,
        "alpha": ALPHA,
        "lambda": LAM,
        "tau_E": TAU_E,
        "late_time_fraction": 0.4,
        "regimes": {
            name: {q: pair(b, q) for q in
                   ("I_pred", "I_mem", "eta", "nu_still", "nu_op", "n_obs")}
            for name, b in blocks.items()
        },
        "I_pred_collapse_ratio_drift_over_stationary": collapse_ratio,
        "eta_range_all_points": {"min": eta_min, "max": eta_max},
        "dpi_violations_I_pred_gt_I_mem": int(dpi_violations),
        "total_sampled_points": int(total_points),
        "tau_w_sweep": {
            "tau_w_values": tau_w_sweep_values,
            "I_pred_mean_late": Ipred_mean_vs_tau_w,
            "I_pred_sem_late": Ipred_sem_vs_tau_w,
            "eta_mean_late": eta_mean_vs_tau_w,
            "eta_sem_late": eta_sem_vs_tau_w,
            "optimal_tau_w_by_I_pred": tau_w_star_Ipred,
            "optimal_tau_w_by_eta": tau_w_star_eta,
        },
        "lemma2_floor": {
            "eta_star_max_over_tau_w": eta_star,
            "c_floor_1_minus_eta_star": c_floor,
            "drift_no_reset_nu_still": nu_still_drift,
            "floor_satisfied_nu_ge_1_minus_c": bool(floor_ok),
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    with SUMMARY_PATH.open("w", encoding="utf-8") as f:
        f.write("=== Markov-drift simulation results (§ 6.1) ===\n")
        f.write("Re-equipped apparatus: eta = I_pred / I_mem,  "
                "I_mem = (K/2) ln(max(N_obs, e))\n")
        f.write("  nu_Still = 1 - eta (memory ballast),  "
                "nu_op = 1 - I_pred/I_opt (secondary)\n")
        f.write(f"k = {K}, K = k(k-1) = {K_PARAMS}, tau_pred = {TAU_PRED}, "
                f"alpha = {ALPHA}, lambda = {LAM}, tau_E = {TAU_E:.0f}\n")
        f.write("All quantities are late-time means over the last 40% of each run "
                "(+- SEM across runs).\n\n")

        order = [
            ("stationary", "Stationary (lam=0, tau_w=inf, T=1e4, 80 runs)"),
            ("drift_no_reset", "Drift, no reset (lam=1e-3, tau_w=inf, T=1e5, 60 runs)"),
            ("drift_reset_100", "Drift + reset tau_w=100 (T=1e5, 60 runs)"),
            ("drift_reset_500", "Drift + reset tau_w=500 (T=1e5, 60 runs)"),
            ("drift_reset_2000", "Drift + reset tau_w=2000 (T=1e5, 60 runs)"),
        ]
        for key, title in order:
            b = blocks[key]
            f.write(f"{title}:\n")
            for q, lbl in (("I_pred", "I_pred"), ("I_mem", "I_mem"), ("eta", "eta"),
                           ("nu_still", "nu_Still"), ("nu_op", "nu_op"),
                           ("n_obs", "N_obs")):
                m, s = b[q]
                f.write(f"    <{lbl:>8}> = {m:.4f} +- {s:.4f}\n")
            f.write("\n")

        f.write(f"I_pred collapse (drift-no-reset / stationary): {collapse_ratio:.4f}\n")
        f.write(f"  (drift {ip_drift:.4f} nats vs stationary {ip_stat:.4f} nats)\n\n")

        f.write(f"eta in [0,1] check over all {total_points} sampled points: "
                f"range [{eta_min:.4f}, {eta_max:.4f}]\n")
        f.write(f"DPI violations (I_pred > I_mem): {dpi_violations}\n\n")

        f.write("Tau_w sweep (lam=1e-3, T=5e4, 40 runs):\n")
        for tw, ipm, ips, etm, ets in zip(
            tau_w_sweep_values, Ipred_mean_vs_tau_w, Ipred_sem_vs_tau_w,
            eta_mean_vs_tau_w, eta_sem_vs_tau_w,
        ):
            f.write(f"  tau_w={tw:>6}: <I_pred>={ipm:.3e}+-{ips:.1e}  "
                    f"<eta>={etm:.3e}+-{ets:.1e}\n")
        f.write(f"  optimal tau_w by <I_pred>: {tau_w_star_Ipred} (tau_E = {TAU_E:.0f})\n")
        f.write(f"  optimal tau_w by <eta>:    {tau_w_star_eta}\n\n")

        f.write("Lemma 2 (memory-ballast floor):\n")
        f.write(f"  eta* = max_tau_w <eta> = {eta_star:.4f}\n")
        f.write(f"  c = 1 - eta* = {c_floor:.4f}\n")
        f.write(f"  drift-no-reset <nu_Still> = {nu_still_drift:.4f} "
                f">= 1 - c = {c_floor:.4f} ? {floor_ok}\n")

    print(f"\nWrote {SUMMARY_PATH} and {SUMMARY_JSON}")


if __name__ == "__main__":
    main()
