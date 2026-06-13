"""ML surrogate «Permuted MNIST» (§ 8.4 main.ru.md, статья #2 Informational Nostalgia).

Честная проверка предсказания Леммы 2 с КОНТРОЛЕМ ВЫУЧЕННОСТИ.
=============================================================
Прежняя метрика nu^op (доля юнитов, чьё зануление не повышает loss на текущей
задаче) оказалась квази-тавтологичной: при быстром дрейфе сеть просто НЕ
ВЫУЧИВАЛА текущую задачу (accuracy ~ случайная), терять было нечего, и nu^op ->1
МЕХАНИЧЕСКИ. То есть nu^op ~ (1 - выученность), а не «доля юнитов, кодирующих
устаревшую старо-permutation структуру». Поэтому здесь:

  1. ACCURACY-GATE / DWELL. dwell = 1/rho_p выбран >> числа шагов до сходимости
     на одной перестановке, так что в основной части скана текущая задача реально
     выучена (accuracy >> случайной = 0.1). Точки с accuracy < ACC_GATE
     помечаются `undertrained` и ИСКЛЮЧАЮТСЯ из вывода о nu (не интерпретируются
     как ностальгия).

  2. ИСТИННАЯ НОСТАЛЬГИЯ nu^op_true. Юнит ностальгичен, если его зануление
     (а) НЕ повышает loss на текущей перестановке выше порога (бесполезен сейчас)
     И (б) ПОВЫШАЕТ loss на хотя бы одной НЕДАВНЕЙ ПРОШЛОЙ перестановке (кодирует
     устаревшую структуру). Это отделяет stale-retained от never-useful/dead.
     Дополнительно отчёт: dead-фракция и useful-фракция.

Запуск:  python main.py     (numpy [+ scikit-learn для load_digits]; matplotlib опционально)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from model import load_dataset, simulate_rho

HERE = Path(__file__).resolve().parent
SUMMARY_TXT = HERE / "results_summary.txt"
SUMMARY_JSON = HERE / "results_summary.json"
LOG_PATH = HERE / "run.log"

# ---------------- Фиксированные параметры ----------------
SEED = 20260524
ETA = 0.1            # learning rate SGD (eta в c = eta/rho_p); крупнее прежнего,
                     # чтобы сходимость на одной перестановке укладывалась в dwell
N_HIDDEN = 32        # скрытых юнитов (единица ablation = один юнит). Берём H=32
                     # (а не 16): при H=16 выученная сеть почти не имеет stale-юнитов
                     # (nu_true ~ 0 во всём learned-диапазоне), сигнал вырожден.
N_STEPS = 120_000    # SGD-шагов на прогон: достаточно, чтобы при самом медленном
                     # rho_p (dwell=1e4) в последней половине было несколько смен
BATCH_SIZE = 32
ABL_THRESHOLD = 0.05  # порог dL для ablation (бесполезен/полезен)
MEASURE_EVERY = 2_000
ABL_SAMPLE = 600     # размер выборки для оценки loss при ablation
N_PAST = 3           # сколько недавних прошлых перестановок держим для теста (б)
N_RUNS = 4           # усреднение по зёрнам -> воспроизводимые числа
ACC_GATE = 0.80      # accuracy-gate: ниже — точка `undertrained`, исключается
                     # из вывода о nu (случайная accuracy для 10 классов = 0.1)

# Сетка rho_p. dwell = 1/rho_p; при eta=0.1 сходимость на одной перестановке
# ~1000 шагов, поэтому при rho_p <= 1e-3 задача выучена (acc >> 0.1).
#   rho_p:  3e-3   1e-3   5e-4   3e-4   1e-4
#   dwell:  333    1000   2000   3333   10000
#   c    :  33     100    200    333    1000
# Точка rho_p=3e-3 (dwell=333) ожидаемо undertrained -> попадёт под gate.
RHO_GRID = [3e-3, 1e-3, 5e-4, 3e-4, 1e-4]


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def main() -> None:
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_lines.append(msg)

    log("== ML surrogate Permuted MNIST (§ 8.4) — честная метрика nu^op_true ==")

    rng_data = np.random.default_rng(SEED)
    X, y, D, C, source = load_dataset(rng_data)
    log(f"датасет: {source}; X = {X.shape}, классов = {C}")
    log(f"eta = {ETA}, H = {N_HIDDEN}, steps = {N_STEPS}, batch = {BATCH_SIZE}, "
        f"abl_threshold = {ABL_THRESHOLD}, n_past = {N_PAST}")
    log(f"N_runs = {N_RUNS}, seed = {SEED}, acc_gate = {ACC_GATE}")
    log(f"rho_p grid = {RHO_GRID}  (c = eta/rho_p)")

    rows = []
    for rho_p in RHO_GRID:
        nus, accs, deads, usefuls, switches = [], [], [], [], []
        per_run = []
        for run in range(N_RUNS):
            rng = np.random.default_rng(SEED + run)
            r = simulate_rho(
                X, y, C,
                rho_p=rho_p, eta=ETA, n_hidden=N_HIDDEN, n_steps=N_STEPS,
                batch_size=BATCH_SIZE, abl_threshold=ABL_THRESHOLD,
                measure_every=MEASURE_EVERY, abl_sample=ABL_SAMPLE,
                n_past=N_PAST, rng=rng,
            )
            nus.append(r["nu_true"])
            accs.append(r["accuracy"])
            deads.append(r["dead_frac"])
            usefuls.append(r["useful_frac"])
            switches.append(r["n_switches"])
            per_run.append({
                "nu_true": r["nu_true"], "accuracy": r["accuracy"],
                "dead_frac": r["dead_frac"], "useful_frac": r["useful_frac"],
                "n_switches": r["n_switches"],
                "n_nu_measurements": r["n_nu_measurements"],
            })
        c_theory = ETA / rho_p
        acc_mean = float(np.nanmean(accs))
        row = {
            "rho_p": rho_p,
            "c_theory": c_theory,
            "dwell": 1.0 / rho_p,
            "accuracy": acc_mean,
            "accuracy_std": float(np.nanstd(accs)),
            "nu_true": float(np.nanmean(nus)),
            "nu_true_std": float(np.nanstd(nus)),
            "dead_frac": float(np.nanmean(deads)),
            "useful_frac": float(np.nanmean(usefuls)),
            "n_switches_mean": float(np.mean(switches)),
            "learned": bool(acc_mean >= ACC_GATE),
            "per_run": per_run,
        }
        rows.append(row)
        flag = "LEARNED " if row["learned"] else "undertr."
        log(f"  rho_p={rho_p:.0e}  c={c_theory:7.1f}  dwell={row['dwell']:7.0f}  "
            f"acc={row['accuracy']:.3f} [{flag}]  "
            f"nu_true={row['nu_true']:.3f}+-{row['nu_true_std']:.3f}  "
            f"dead={row['dead_frac']:.3f}  useful={row['useful_frac']:.3f}")

    # --- проверка предсказания ТОЛЬКО на выученных точках (accuracy-gate) ---
    learned_rows = [r for r in rows if r["learned"]]
    log("")
    log("== Проверка предсказания (accuracy-gated) ==")
    log(f"выученных точек (acc >= {ACC_GATE}): {len(learned_rows)} из {len(rows)}")

    prediction_check: dict = {"acc_gate": ACC_GATE,
                              "n_learned_points": len(learned_rows)}

    if len(learned_rows) >= 2:
        rho_arr = np.array([r["rho_p"] for r in learned_rows])
        nu_arr = np.array([r["nu_true"] for r in learned_rows])
        acc_arr = np.array([r["accuracy"] for r in learned_rows])
        c_arr = np.array([r["c_theory"] for r in learned_rows])

        sp_rho_nu = spearman(rho_arr, nu_arr)
        sp_acc_nu = spearman(acc_arr, nu_arr)
        sp_c_nu = spearman(c_arr, nu_arr)
        # Лемма 2: ностальгия растёт с частотой дрейфа rho_p (Спирмен ~ +1)
        log(f"Спирмен(rho_p, nu_true)   = {sp_rho_nu:+.3f}  "
            f"(Лемма 2 ожидает > 0: ностальгия растёт с частотой дрейфа)")
        log(f"Спирмен(c,     nu_true)   = {sp_c_nu:+.3f}  "
            f"(зеркально ожидается < 0)")
        # конфаунд: тривиальная антикорреляция nu с accuracy ДОЛЖНА исчезнуть
        log(f"Спирмен(accuracy, nu_true)= {sp_acc_nu:+.3f}  "
            f"(прежний конфаунд: было ~ -1; теперь не должно быть тривиальной "
            f"антикорреляции)")
        prediction_check.update({
            "spearman_rho_nu_learned": sp_rho_nu,
            "spearman_c_nu_learned": sp_c_nu,
            "spearman_acc_nu_learned": sp_acc_nu,
        })
    else:
        log("Недостаточно выученных точек для корреляционного вывода.")

    summary = {
        "experiment": "ML surrogate Permuted MNIST (article #2 § 8.4), nu^op_true",
        "dataset_source": source,
        "X_shape": list(X.shape),
        "n_classes": C,
        "eta": ETA,
        "n_hidden": N_HIDDEN,
        "n_steps": N_STEPS,
        "batch_size": BATCH_SIZE,
        "abl_threshold": ABL_THRESHOLD,
        "measure_every": MEASURE_EVERY,
        "abl_sample": ABL_SAMPLE,
        "n_past": N_PAST,
        "n_runs": N_RUNS,
        "seed": SEED,
        "acc_gate": ACC_GATE,
        "rho_grid": RHO_GRID,
        "scan": rows,
        "prediction_check": prediction_check,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    with SUMMARY_TXT.open("w", encoding="utf-8") as f:
        f.write("=== ML surrogate Permuted MNIST (статья #2, § 8.4) — nu^op_true ===\n")
        f.write(f"датасет: {source}\n")
        f.write(f"X = {tuple(X.shape)}, классов = {C}\n")
        f.write(f"eta = {ETA}, H = {N_HIDDEN}, steps = {N_STEPS}, batch = {BATCH_SIZE}\n")
        f.write(f"abl_threshold = {ABL_THRESHOLD}, abl_sample = {ABL_SAMPLE}, "
                f"measure_every = {MEASURE_EVERY}, n_past = {N_PAST}\n")
        f.write(f"N_runs = {N_RUNS}, seed = {SEED}, acc_gate = {ACC_GATE}\n\n")

        f.write("--- Скан rho_p (c = eta/rho_p), КОНТРОЛЬ ВЫУЧЕННОСТИ ---\n")
        f.write(f"{'rho_p':>8s}  {'c':>7s}  {'dwell':>7s}  {'acc':>6s}  "
                f"{'learned':>7s}  {'nu_true':>8s}  {'dead':>6s}  {'useful':>6s}\n")
        for r in rows:
            f.write(f"{r['rho_p']:>8.1e}  {r['c_theory']:>7.1f}  "
                    f"{r['dwell']:>7.0f}  {r['accuracy']:>6.3f}  "
                    f"{('yes' if r['learned'] else 'NO'):>7s}  "
                    f"{r['nu_true']:>8.3f}  {r['dead_frac']:>6.3f}  "
                    f"{r['useful_frac']:>6.3f}\n")

        f.write("\n--- Проверка предсказания (только выученные точки) ---\n")
        f.write(f"выученных точек (acc >= {ACC_GATE}): "
                f"{len(learned_rows)} из {len(rows)}\n")
        if "spearman_rho_nu_learned" in prediction_check:
            f.write(f"Спирмен(rho_p, nu_true)    = "
                    f"{prediction_check['spearman_rho_nu_learned']:+.3f}  "
                    f"(Лемма 2 ожидает > 0)\n")
            f.write(f"Спирмен(c,     nu_true)    = "
                    f"{prediction_check['spearman_c_nu_learned']:+.3f}\n")
            f.write(f"Спирмен(accuracy, nu_true) = "
                    f"{prediction_check['spearman_acc_nu_learned']:+.3f}  "
                    f"(конфаунд с выученностью должен исчезнуть)\n")
        f.write("\nИнтерпретация — см. README.ru.md / README.md «Честная рамка».\n")

    LOG_PATH.write_text("\n".join(log_lines), encoding="utf-8")

    # --- опциональная фигура ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        c_all = np.array([r["c_theory"] for r in rows])
        nu_all = np.array([r["nu_true"] for r in rows])
        acc_all = np.array([r["accuracy"] for r in rows])
        learned_mask = np.array([r["learned"] for r in rows])

        fig, ax = plt.subplots(figsize=(7.8, 5.0))
        ax.errorbar(c_all[learned_mask], nu_all[learned_mask],
                    yerr=np.array([r["nu_true_std"] for r in rows])[learned_mask],
                    fmt="o-", color="tab:blue", ms=8, lw=1.5, capsize=4,
                    label=r"$\nu^{\rm op}_{\rm true}$ (learned, acc-gated)")
        if (~learned_mask).any():
            ax.scatter(c_all[~learned_mask], nu_all[~learned_mask],
                       marker="x", color="tab:red", s=60,
                       label="undertrained (excluded)")
        ax2 = ax.twinx()
        ax2.plot(c_all, acc_all, "s--", color="tab:green", ms=5, lw=1.0,
                 alpha=0.6, label="accuracy")
        ax2.axhline(ACC_GATE, color="tab:green", ls=":", lw=0.8, alpha=0.6)
        ax2.set_ylabel("current-task accuracy", color="tab:green")
        ax.set_xscale("log")
        ax.set_xlabel(r"$c = \eta/\rho_p$  (slower switching $\to$ larger $c$)")
        ax.set_ylabel(r"true nostalgic fraction $\nu^{\rm op}_{\rm true}$",
                      color="tab:blue")
        ax.set_title(r"Permuted MNIST: $\nu^{\rm op}_{\rm true}$ vs "
                     r"$c=\eta/\rho_p$ under learnedness control")
        for r in rows:
            ax.annotate(rf"$\rho_p={r['rho_p']:.0e}$", (r["c_theory"], r["nu_true"]),
                        textcoords="offset points", xytext=(5, 6), fontsize=8)
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc="best", fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        fig_path = HERE / "fig_nu_vs_c.png"
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        log(f"saved {fig_path}")
        LOG_PATH.write_text("\n".join(log_lines), encoding="utf-8")
    except Exception as exc:  # matplotlib опционален
        log(f"(matplotlib недоступен, фигура пропущена: {exc})")

    print(f"\nWrote {SUMMARY_TXT}, {SUMMARY_JSON}, {LOG_PATH}")


if __name__ == "__main__":
    main()
