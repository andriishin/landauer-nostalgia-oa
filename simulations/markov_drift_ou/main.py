"""Run OU-class simulation for paper § 6.2 and produce figs 4, 5.

Re-equipped apparatus: efficiency eta = I_pred / I_mem (was I_pred / N_max).
The old "slope K/(2R) = 28" conjecture fit was an artefact of the energy-budget
denominator N_max = R*t and has been REMOVED.

The central result of § 6.2 is a CONTRAST between two online learners run on the
SAME drifting trajectory:
  - Robbins-Monro (decaying step c/sqrt(t)): FREEZES under drift -- its step
    vanishes, the estimate stops tracking, and the operational nostalgia nu_op
    climbs toward 1.
  - constant step (s = 0.3, validated optimal at sigma = 0.1): TRACKS the moving
    target indefinitely, paying a residual lag/noise cost so nu_op settles to a
    finite FLOOR well below the frozen RM value.
This is a design lesson, not a no-go theorem: under non-stationarity the optimal
estimator is not the asymptotically-optimal stationary one.

Usage:
    python main.py

Outputs (descriptive PNG preview + vector PDF for the manuscript/README):
    - paper/figs/fig4_ou_nuop_vs_t.png            (+ Fig4.pdf)
        Fig 4 (carrier): nu_op(t) of BOTH learners under drift -- RM freezing
        (rising toward ~0.86) vs constant step holding a floor (~0.34); inset
        shows I_pred(t) of both.
    - paper/figs/fig5_ou_optimal_step.png   (+ Fig5.pdf)
        Fig 5 (design decides): constant-step nu_op vs step magnitude at
        sigma = 0.1 -- an interior optimum (too small -> lag, too large ->
        noise), with the frozen-RM level as a horizontal reference line.
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

# ----- Global parameters (see § 5.2 / § 6.2 of main.ru.md) -----
# Note on sigma choice: with very small sigma the OU stationary variance
# sigma^2/(2 lambda) keeps theta near zero, so the transition matrix is
# indistinguishable from uniform and I_opt collapses to noise. We use
# sigma = 0.1 (OU stationary variance ~ 5) as the central working point, which
# yields a genuinely informative chain (I_opt ~ 0.7 nats per step). A sweep over
# sigma in {0.03, 0.06, 0.1, 0.2, 0.4} probes how the drift strength sets the
# floor of the operational nostalgia nu_op.
K_STATES = 8                # number of states
K_PARAMS = K_STATES * (K_STATES - 1)   # = 56 off-diagonal logits
LAM = 1e-3                  # OU mean-reversion rate
SIGMA = 0.1                 # OU diffusion amplitude (central working point)
T_TOTAL = 20_000            # simulation length in steps
N_RUNS = 50                 # independent OU realisations
THETA0_STD = 0.1            # std of initial logits
TAU_PRED = 1                # one-step predictive horizon
MEASURE_EVERY = 200         # subsampling for measurement
LEARN_C = 3.0               # initial learning rate of online ML (Robbins-Monro c)
CONST_STEP = 0.3            # constant-step learner working point (optimal at sigma=0.1)
SEED = 20260524

# sigma sweep (drift-strength dependence of nu_op floor for both learners)
SIGMA_SWEEP = [0.03, 0.06, 0.1, 0.2, 0.4]
SWEEP_RUNS = 40

# constant-step magnitude sweep for fig 5 (interior optimum at central sigma)
STEP_SWEEP = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0]
STEP_SWEEP_RUNS = 40


def late_time_mean(runs: list[dict], key: str, frac: float = 0.5) -> tuple[float, float]:
    """Per-run late-time average of `key` over the last `frac` of the run, then
    mean +- SEM across runs. Returns (mean, sem)."""
    T = float(runs[0]["times"][-1])
    per_run = []
    for r in runs:
        m = r["times"] >= T * (1.0 - frac)
        if m.any():
            per_run.append(float(np.mean(r[key][m])))
    if not per_run:
        return float("nan"), float("nan")
    arr = np.array(per_run)
    return float(arr.mean()), float(arr.std() / np.sqrt(len(arr)))


def run_condition(lam: float, sigma: float, T: int, n_runs: int,
                  measure_every: int, seed_offset: int,
                  const_step: float = CONST_STEP) -> dict:
    runs = []
    for r in range(n_runs):
        rng = np.random.default_rng(SEED + seed_offset + r)
        run = simulate_run(
            k=K_STATES,
            T=T,
            lam=lam,
            sigma=sigma,
            tau_pred=TAU_PRED,
            rng=rng,
            measure_every=measure_every,
            theta0_std=THETA0_STD,
            learn_c=LEARN_C,
            const_step=const_step,
        )
        runs.append(run)
    avg = average_runs(runs)
    avg["raw"] = runs
    return avg


def main() -> None:
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_lines.append(msg)

    log("== OU-drift § 6.2: TWO learners (Robbins-Monro vs constant step) ==")
    log(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}")
    log(f"lambda = {LAM}, sigma = {SIGMA}, T = {T_TOTAL}, n_runs = {N_RUNS}")
    log(f"theta0_std = {THETA0_STD}, tau_pred = {TAU_PRED}, measure_every = {MEASURE_EVERY}")
    log(f"Robbins-Monro step = c/sqrt(t), c = {LEARN_C};  constant step s = {CONST_STEP}")
    log(f"OU stationary variance sigma^2/(2 lambda) = {SIGMA**2 / (2 * LAM):.3e}")
    log(f"Adiabatic check: lambda * tau_relax ~ lambda * O(1) = {LAM:.1e} (<< 1: OK)")
    log("I_mem(t) = (K/2) ln(max(N_obs, e)),  N_obs = t (online learner, no window)")
    log("eta = I_pred/I_mem ;  nu_op = 1 - I_pred/I_opt (canonical)")

    t0 = time.time()

    # ---------- Central working point (both learners, same trajectory) ----------
    avg = run_condition(lam=LAM, sigma=SIGMA, T=T_TOTAL, n_runs=N_RUNS,
                        measure_every=MEASURE_EVERY, seed_offset=0)
    runs = avg["raw"]
    times = avg["times"]
    log(f"Central working point done (elapsed {time.time() - t0:.1f}s)")

    # ---------- sigma sweep (both learners) ----------
    log("\n== sigma sweep (drift strength): RM flat-high vs const responsive ==")
    sw_nuop_rm_m, sw_nuop_rm_s = [], []
    sw_nuop_c_m, sw_nuop_c_s = [], []
    sw_ipred_rm_m, sw_ipred_c_m = [], []
    sweep_avgs = {}
    for i, sig in enumerate(SIGMA_SWEEP):
        cond = run_condition(lam=LAM, sigma=sig, T=T_TOTAL, n_runs=SWEEP_RUNS,
                            measure_every=MEASURE_EVERY, seed_offset=1000 + 1000 * i)
        sweep_avgs[sig] = cond
        no_rm_m, no_rm_s = late_time_mean(cond["raw"], "nu_op_rm")
        no_c_m, no_c_s = late_time_mean(cond["raw"], "nu_op_const")
        ip_rm_m, _ = late_time_mean(cond["raw"], "I_pred_rm")
        ip_c_m, _ = late_time_mean(cond["raw"], "I_pred_const")
        sw_nuop_rm_m.append(no_rm_m); sw_nuop_rm_s.append(no_rm_s)
        sw_nuop_c_m.append(no_c_m); sw_nuop_c_s.append(no_c_s)
        sw_ipred_rm_m.append(ip_rm_m); sw_ipred_c_m.append(ip_c_m)
        log(f"  sigma={sig:>5}: RM <nu_op>={no_rm_m:.3f}+-{no_rm_s:.3f}  "
            f"const <nu_op>={no_c_m:.3f}+-{no_c_s:.3f}")

    # ---------- constant-step magnitude sweep (Fig 5: interior optimum) ----------
    log(f"\n== constant-step magnitude sweep at sigma={SIGMA} (Fig 5) ==")
    step_nuop_m, step_nuop_s = [], []
    step_ipred_m = []
    for i, s in enumerate(STEP_SWEEP):
        cond = run_condition(lam=LAM, sigma=SIGMA, T=T_TOTAL, n_runs=STEP_SWEEP_RUNS,
                            measure_every=MEASURE_EVERY, seed_offset=5000 + 1000 * i,
                            const_step=s)
        no_m, no_s = late_time_mean(cond["raw"], "nu_op_const")
        ip_m, _ = late_time_mean(cond["raw"], "I_pred_const")
        step_nuop_m.append(no_m); step_nuop_s.append(no_s)
        step_ipred_m.append(ip_m)
        log(f"  step={s:>5}: const <nu_op>={no_m:.3f}+-{no_s:.3f}  "
            f"<I_pred>={ip_m:.3f}")
    step_arr = np.array(STEP_SWEEP, dtype=float)
    best_idx = int(np.argmin(step_nuop_m))
    best_step = STEP_SWEEP[best_idx]
    best_step_nuop = step_nuop_m[best_idx]
    log(f"  optimal constant step ~ {best_step} (min <nu_op>={best_step_nuop:.3f})")

    log(f"\nTotal simulation time: {time.time() - t0:.1f}s")

    # ---------- Late-time diagnostics (central working point), per learner ----------
    Iopt_m, Iopt_s = late_time_mean(runs, "I_opt")
    Imem_m, Imem_s = late_time_mean(runs, "I_mem")
    nobs_m, _ = late_time_mean(runs, "n_obs")

    # Robbins-Monro
    Ipred_rm_m, Ipred_rm_s = late_time_mean(runs, "I_pred_rm")
    eta_rm_m, eta_rm_s = late_time_mean(runs, "eta_rm")
    nustill_rm_m, nustill_rm_s = late_time_mean(runs, "nu_still_rm")
    nuop_rm_m, nuop_rm_s = late_time_mean(runs, "nu_op_rm")
    thetaerr_rm_m, thetaerr_rm_s = late_time_mean(runs, "theta_err_rm")
    # constant step
    Ipred_c_m, Ipred_c_s = late_time_mean(runs, "I_pred_const")
    eta_c_m, eta_c_s = late_time_mean(runs, "eta_const")
    nustill_c_m, nustill_c_s = late_time_mean(runs, "nu_still_const")
    nuop_c_m, nuop_c_s = late_time_mean(runs, "nu_op_const")
    thetaerr_c_m, thetaerr_c_s = late_time_mean(runs, "theta_err_const")

    log("\n== Late-time diagnostics (central, last 50% of run) ==")
    log(f"  <I_opt> = {Iopt_m:.4f}+-{Iopt_s:.4f}  <I_mem> = {Imem_m:.4f}+-{Imem_s:.4f}  "
        f"<N_obs> = {nobs_m:.0f}")
    log("  -- Robbins-Monro (decaying step c/sqrt(t)) --")
    log(f"    <I_pred>  = {Ipred_rm_m:.4f} +- {Ipred_rm_s:.4f} nats")
    log(f"    <eta>     = {eta_rm_m:.4f} +- {eta_rm_s:.4f}")
    log(f"    <nu_op>   = {nuop_rm_m:.4f} +- {nuop_rm_s:.4f}   (CANONICAL; should be HIGH/frozen)")
    log(f"    <theta_err> = {thetaerr_rm_m:.4f} +- {thetaerr_rm_s:.4f}")
    log(f"  -- constant step s = {CONST_STEP} --")
    log(f"    <I_pred>  = {Ipred_c_m:.4f} +- {Ipred_c_s:.4f} nats")
    log(f"    <eta>     = {eta_c_m:.4f} +- {eta_c_s:.4f}")
    log(f"    <nu_op>   = {nuop_c_m:.4f} +- {nuop_c_s:.4f}   (CANONICAL; should be a low FLOOR)")
    log(f"    <theta_err> = {thetaerr_c_m:.4f} +- {thetaerr_c_s:.4f}")

    # ---------- RM freezing demonstration (early / mid / late nu_op) ----------
    n = len(times)
    i_early = max(1, n // 10)        # ~10% in
    i_mid = n // 2
    i_late = n - 1
    nuop_rm_curve = avg["nu_op_rm"]
    nuop_c_curve = avg["nu_op_const"]
    log("\n== RM freezing vs const tracking (nu_op early / mid / late) ==")
    log(f"  t = {times[i_early]:.0f}: RM nu_op = {nuop_rm_curve[i_early]:.3f}, "
        f"const = {nuop_c_curve[i_early]:.3f}")
    log(f"  t = {times[i_mid]:.0f}: RM nu_op = {nuop_rm_curve[i_mid]:.3f}, "
        f"const = {nuop_c_curve[i_mid]:.3f}")
    log(f"  t = {times[i_late]:.0f}: RM nu_op = {nuop_rm_curve[i_late]:.3f}, "
        f"const = {nuop_c_curve[i_late]:.3f}")
    log(f"  => RM nu_op rises {nuop_rm_curve[i_early]:.3f} -> {nuop_rm_curve[i_late]:.3f} "
        f"(FREEZING); const stays ~{nuop_c_curve[i_late]:.3f} (tracking floor)")

    # ---------- Bound checks over ALL sampled points (BOTH learners) ----------
    all_runs = list(runs)
    for sig in SIGMA_SWEEP:
        all_runs += sweep_avgs[sig]["raw"]
    eta_min, eta_max = np.inf, -np.inf
    nuop_min, nuop_max = np.inf, -np.inf
    dpi_violations = 0       # I_pred > I_mem  (eta > 1)  -- either learner
    iopt_violations = 0      # I_pred > I_opt  (nu_op < 0) -- either learner
    total_points = 0
    for r in all_runs:
        for tag in ("rm", "const"):
            e = r[f"eta_{tag}"]
            no = r[f"nu_op_{tag}"]
            eta_min = min(eta_min, float(e.min())); eta_max = max(eta_max, float(e.max()))
            nuop_min = min(nuop_min, float(no.min())); nuop_max = max(nuop_max, float(no.max()))
            dpi_violations += int(np.sum(r[f"I_pred_{tag}"] > r["I_mem"] + 1e-9))
            iopt_violations += int(np.sum(r[f"I_pred_{tag}"] > r["I_opt"] + 1e-9))
            total_points += e.size

    log("\n== Bound checks over ALL sampled points (both learners) ==")
    log(f"  total sampled points: {total_points}")
    log(f"  eta range:   [{eta_min:.4f}, {eta_max:.4f}]  "
        f"(DPI violations I_pred>I_mem: {dpi_violations})")
    log(f"  nu_op range: [{nuop_min:.4f}, {nuop_max:.4f}]  "
        f"(I_pred>I_opt violations: {iopt_violations})")

    # ---------- FIG 4 (carrier): nu_op(t) of BOTH learners under drift ----------
    # RM freezes (nu_op climbs); constant step holds a floor. Inset: I_pred(t).
    nuop_rm = avg["nu_op_rm"]; nuop_rm_std = avg["nu_op_rm_std"]
    nuop_c = avg["nu_op_const"]; nuop_c_std = avg["nu_op_const_std"]
    ipred_rm = avg["I_pred_rm"]
    ipred_c = avg["I_pred_const"]
    I_opt_curve = avg["I_opt"]

    fig, ax = plt.subplots(figsize=(8.2, 5.0))

    ax.plot(times, nuop_rm, color="tab:red", lw=1.9,
            label=r"Robbins--Monro $c/\sqrt{t}$ (freezes)")
    ax.fill_between(times, np.clip(nuop_rm - nuop_rm_std, 0, 1),
                    np.clip(nuop_rm + nuop_rm_std, 0, 1),
                    alpha=0.15, color="tab:red")
    ax.plot(times, nuop_c, color="tab:blue", lw=1.9,
            label=fr"constant step $s={CONST_STEP:g}$ (tracks)")
    ax.fill_between(times, np.clip(nuop_c - nuop_c_std, 0, 1),
                    np.clip(nuop_c + nuop_c_std, 0, 1),
                    alpha=0.15, color="tab:blue")
    ax.axhline(nuop_rm_m, ls=":", color="tab:red", alpha=0.55,
               label=fr"RM late-time $\langle\nu^{{\mathrm{{op}}}}\rangle={nuop_rm_m:.2f}$ (frozen)")
    ax.axhline(nuop_c_m, ls=":", color="tab:blue", alpha=0.55,
               label=fr"const floor $\langle\nu^{{\mathrm{{op}}}}\rangle={nuop_c_m:.2f}$")
    ax.set_xlabel("time $t$ (steps)")
    ax.set_ylabel(r"operational nostalgia $\nu^{\mathrm{op}}=1-I_{\mathrm{pred}}/I_{\mathrm{opt}}$")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_title(r"Freezing vs tracking under drift: $\nu^{\mathrm{op}}(t)$ "
                 rf"($k={K_STATES}$, $\lambda={LAM:g}$, $\sigma={SIGMA:g}$)")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # inset: I_pred(t) of both, with I_opt reference.
    # Placed upper-left over the (shaded) confidence band so it clears the
    # constant-step nu_op line (~0.3, full x) and the RM solid line (<0.56 on
    # the left half) -- avoids hiding either carrier curve.
    axin = ax.inset_axes([0.14, 0.57, 0.40, 0.40])
    axin.plot(times, ipred_rm, color="tab:red", lw=1.3,
              label="RM")
    axin.plot(times, ipred_c, color="tab:blue", lw=1.3,
              label="const")
    axin.plot(times, I_opt_curve, color="tab:grey", lw=1.0, ls=":",
              label=r"$I_{\mathrm{opt}}$")
    axin.set_xscale("log")
    axin.set_xlabel(r"$t$", fontsize=7)
    axin.set_ylabel(r"$I_{\mathrm{pred}}$ (nats)", fontsize=7)
    axin.tick_params(labelsize=6)
    axin.legend(fontsize=6, loc="upper right")
    axin.grid(True, which="both", alpha=0.25)

    fig.tight_layout()
    fig4_path = PAPER_FIGS / "fig4_ou_nuop_vs_t.png"  # descriptive name matches content
    fig.savefig(fig4_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig4.pdf")
    plt.close(fig)
    log(f"saved {fig4_path}")

    # ---------- FIG 5 (design decides): const-step nu_op vs step magnitude -----
    # Interior optimum (small step -> lag, large step -> noise) beats frozen RM.
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    nu_step = np.array(step_nuop_m); nu_step_err = np.array(step_nuop_s)

    ax.errorbar(step_arr, nu_step, yerr=nu_step_err, marker="o", capsize=3, lw=1.7,
                color="tab:blue",
                label=r"constant step $\langle\nu^{\mathrm{op}}\rangle$ (floor)")
    ax.axhline(nuop_rm_m, ls="--", color="tab:red", lw=1.5,
               label=fr"frozen Robbins--Monro $\langle\nu^{{\mathrm{{op}}}}\rangle={nuop_rm_m:.2f}$")
    # mark the interior optimum
    ax.scatter([best_step], [best_step_nuop], s=90, facecolors="none",
               edgecolors="tab:green", lw=2.0, zorder=5,
               label=fr"optimal step $s^\star\approx{best_step:g}$ "
                     fr"($\langle\nu^{{\mathrm{{op}}}}\rangle={best_step_nuop:.2f}$)")
    ax.set_xscale("log")
    ax.set_xlabel(r"constant step size $s$")
    ax.set_ylabel(r"late-time $\langle\nu^{\mathrm{op}}\rangle = "
                  r"\langle 1 - I_{\mathrm{pred}}/I_{\mathrm{opt}}\rangle$")
    ax.set_ylim(0, 1.02)
    ax.set_title(r"Design decides: interior optimum in step size beats frozen RM "
                 rf"($\sigma={SIGMA:g}$, $\lambda={LAM:g}$, $k={K_STATES}$)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig5_path = PAPER_FIGS / "fig5_ou_optimal_step.png"  # descriptive name matches content
    fig.savefig(fig5_path, dpi=150)
    fig.savefig(PAPER_FIGS / "Fig5.pdf")
    plt.close(fig)
    log(f"saved {fig5_path}")

    # ---------- Write summary ----------
    summary = {
        "apparatus": "TWO learners on one trajectory. eta = I_pred / I_mem ;  "
                     "I_mem = (K/2) ln(max(N_obs,e)), N_obs=t ;  "
                     "nu_Still = 1 - eta (secondary) ;  "
                     "nu_op = 1 - I_pred/I_opt (canonical nostalgia). "
                     "rm = Robbins-Monro c/sqrt(t) (freezes) ;  "
                     "const = constant step (tracks to a floor).",
        "k": K_STATES,
        "K_params": K_PARAMS,
        "lambda": LAM,
        "sigma": SIGMA,
        "T_total": T_TOTAL,
        "n_runs": N_RUNS,
        "tau_pred": TAU_PRED,
        "theta0_std": THETA0_STD,
        "learn_c": LEARN_C,
        "const_step": CONST_STEP,
        "measure_every": MEASURE_EVERY,
        "seed": SEED,
        "late_time_fraction": 0.5,
        "central_late_time": {
            "shared": {
                "I_opt": {"mean": Iopt_m, "sem": Iopt_s},
                "I_mem": {"mean": Imem_m, "sem": Imem_s},
                "n_obs": nobs_m,
            },
            "robbins_monro": {
                "I_pred": {"mean": Ipred_rm_m, "sem": Ipred_rm_s},
                "eta": {"mean": eta_rm_m, "sem": eta_rm_s},
                "nu_still": {"mean": nustill_rm_m, "sem": nustill_rm_s},
                "nu_op": {"mean": nuop_rm_m, "sem": nuop_rm_s},
                "theta_err": {"mean": thetaerr_rm_m, "sem": thetaerr_rm_s},
            },
            "constant_step": {
                "step": CONST_STEP,
                "I_pred": {"mean": Ipred_c_m, "sem": Ipred_c_s},
                "eta": {"mean": eta_c_m, "sem": eta_c_s},
                "nu_still": {"mean": nustill_c_m, "sem": nustill_c_s},
                "nu_op": {"mean": nuop_c_m, "sem": nuop_c_s},
                "theta_err": {"mean": thetaerr_c_m, "sem": thetaerr_c_s},
            },
        },
        "rm_freezing": {
            "times": [float(times[i_early]), float(times[i_mid]), float(times[i_late])],
            "nu_op_rm": [float(nuop_rm_curve[i_early]), float(nuop_rm_curve[i_mid]),
                         float(nuop_rm_curve[i_late])],
            "nu_op_const": [float(nuop_c_curve[i_early]), float(nuop_c_curve[i_mid]),
                            float(nuop_c_curve[i_late])],
        },
        "bound_checks": {
            "total_sampled_points": int(total_points),
            "eta_range": {"min": float(eta_min), "max": float(eta_max)},
            "dpi_violations_I_pred_gt_I_mem": int(dpi_violations),
            "nu_op_range": {"min": float(nuop_min), "max": float(nuop_max)},
            "iopt_violations_I_pred_gt_I_opt": int(iopt_violations),
        },
        "sigma_sweep": {
            "sigma_values": SIGMA_SWEEP,
            "n_runs_per_sigma": SWEEP_RUNS,
            "nu_op_rm_mean_late": sw_nuop_rm_m,
            "nu_op_rm_sem_late": sw_nuop_rm_s,
            "nu_op_const_mean_late": sw_nuop_c_m,
            "nu_op_const_sem_late": sw_nuop_c_s,
            "I_pred_rm_mean_late": sw_ipred_rm_m,
            "I_pred_const_mean_late": sw_ipred_c_m,
        },
        "step_sweep": {
            "step_values": STEP_SWEEP,
            "n_runs_per_step": STEP_SWEEP_RUNS,
            "sigma": SIGMA,
            "nu_op_const_mean_late": step_nuop_m,
            "nu_op_const_sem_late": step_nuop_s,
            "I_pred_const_mean_late": step_ipred_m,
            "optimal_step": best_step,
            "optimal_step_nu_op": best_step_nuop,
            "rm_nu_op_reference": nuop_rm_m,
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    with SUMMARY_PATH.open("w", encoding="utf-8") as f:
        f.write("=== OU-drift simulation results (§ 6.2): TWO learners ===\n")
        f.write("Apparatus: eta = I_pred / I_mem,  I_mem = (K/2) ln(max(N_obs, e)), N_obs = t\n")
        f.write("  nu_op = 1 - I_pred/I_opt (CANONICAL nostalgia),  nu_Still = 1 - eta (secondary)\n")
        f.write("  rm   = Robbins-Monro decaying step c/sqrt(t)  -> FREEZES under drift\n")
        f.write("  const= constant step s (validated optimal at sigma=0.1) -> TRACKS to a floor\n")
        f.write(f"k = {K_STATES}, K = k(k-1) = {K_PARAMS}, tau_pred = {TAU_PRED}\n")
        f.write(f"lambda = {LAM}, sigma = {SIGMA}, T = {T_TOTAL}, n_runs = {N_RUNS}\n")
        f.write(f"theta0_std = {THETA0_STD}, learn_c = {LEARN_C}, const_step = {CONST_STEP}, "
                f"seed = {SEED}\n")
        f.write(f"OU stationary variance sigma^2/(2 lambda) = {SIGMA**2 / (2 * LAM):.3e}\n\n")

        f.write("--- Central working point, late-time means (last 50%, +- SEM across runs) ---\n")
        f.write(f"  shared:   <I_opt> = {Iopt_m:.4f}+-{Iopt_s:.4f}   "
                f"<I_mem> = {Imem_m:.4f}+-{Imem_s:.4f}   <N_obs> = {nobs_m:.0f}\n")
        f.write("  Robbins-Monro (decaying step c/sqrt(t)):\n")
        f.write(f"    <I_pred>  = {Ipred_rm_m:.4f} +- {Ipred_rm_s:.4f} nats\n")
        f.write(f"    <eta>     = {eta_rm_m:.4f} +- {eta_rm_s:.4f}\n")
        f.write(f"    <nu_op>   = {nuop_rm_m:.4f} +- {nuop_rm_s:.4f}  (CANONICAL; HIGH = frozen)\n")
        f.write(f"    <theta_err> = {thetaerr_rm_m:.4f} +- {thetaerr_rm_s:.4f}\n")
        f.write(f"  constant step s = {CONST_STEP}:\n")
        f.write(f"    <I_pred>  = {Ipred_c_m:.4f} +- {Ipred_c_s:.4f} nats\n")
        f.write(f"    <eta>     = {eta_c_m:.4f} +- {eta_c_s:.4f}\n")
        f.write(f"    <nu_op>   = {nuop_c_m:.4f} +- {nuop_c_s:.4f}  (CANONICAL; LOW = tracking floor)\n")
        f.write(f"    <theta_err> = {thetaerr_c_m:.4f} +- {thetaerr_c_s:.4f}\n\n")

        f.write("--- RM freezing vs const tracking (nu_op early / mid / late) ---\n")
        f.write(f"  t={times[i_early]:.0f}: RM={nuop_rm_curve[i_early]:.3f}  "
                f"const={nuop_c_curve[i_early]:.3f}\n")
        f.write(f"  t={times[i_mid]:.0f}: RM={nuop_rm_curve[i_mid]:.3f}  "
                f"const={nuop_c_curve[i_mid]:.3f}\n")
        f.write(f"  t={times[i_late]:.0f}: RM={nuop_rm_curve[i_late]:.3f}  "
                f"const={nuop_c_curve[i_late]:.3f}\n")
        f.write(f"  => RM nu_op rises {nuop_rm_curve[i_early]:.3f} -> {nuop_rm_curve[i_late]:.3f} "
                f"(FREEZING); const stays ~{nuop_c_curve[i_late]:.3f} (floor)\n\n")

        f.write("--- Bound checks over ALL sampled points (both learners) ---\n")
        f.write(f"  total sampled points: {total_points}\n")
        f.write(f"  eta in [0,1]:    range [{eta_min:.4f}, {eta_max:.4f}]  "
                f"(I_pred>I_mem violations: {dpi_violations})\n")
        f.write(f"  nu_op in [0,1]:  range [{nuop_min:.4f}, {nuop_max:.4f}]  "
                f"(I_pred>I_opt violations: {iopt_violations})\n\n")

        f.write("--- sigma sweep (drift strength), late-time <nu_op> ---\n")
        f.write("  (RM stays flat-high ~ frozen; const responds to drift)\n")
        for sig, no_rm, no_c in zip(SIGMA_SWEEP, sw_nuop_rm_m, sw_nuop_c_m):
            f.write(f"  sigma={sig:>5}: RM <nu_op>={no_rm:.3f}   const <nu_op>={no_c:.3f}\n")
        f.write("\n")

        f.write(f"--- constant-step magnitude sweep at sigma={SIGMA} (Fig 5: interior optimum) ---\n")
        for s, no_m, no_s in zip(STEP_SWEEP, step_nuop_m, step_nuop_s):
            mark = "  <-- optimum" if s == best_step else ""
            f.write(f"  step={s:>5}: const <nu_op>={no_m:.3f}+-{no_s:.3f}{mark}\n")
        f.write(f"  optimal step ~ {best_step} (<nu_op>={best_step_nuop:.3f}) "
                f"vs frozen RM <nu_op>={nuop_rm_m:.3f}\n")

    LOG_PATH.write_text("\n".join(log_lines), encoding="utf-8")
    print(f"\nWrote {SUMMARY_PATH}, {SUMMARY_JSON}, {LOG_PATH}")


if __name__ == "__main__":
    main()
