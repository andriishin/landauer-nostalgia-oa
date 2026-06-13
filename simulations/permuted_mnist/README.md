# ML paradigm case: Permuted MNIST (§ 8.4)

## What it demonstrates

A paradigm-case ML surrogate for the falsifiable prediction of § 8.4 of
`paper/main.ru.md` (the ML instantiation of Lemma 2), tested **with an explicit
learnedness control**. The prediction under test: the *true* nostalgic fraction
`nu^op_true` — units that retain stale, once-predictive structure — should grow
with the permutation-switch rate `rho_p` (smaller canonical `c = eta/rho_p`).

**This surrogate reports an honest negative/reversed result on that naive
reading** (see *Honest framing*): once learnedness is controlled, the monotone
`nu` ↑ with `rho_p` that a naive ablation metric would show turns out to be an
artifact, and the surrogate instead exhibits a different, well-defined effect.

The canonical parameter being scanned:

```
c := M_refresh_dot * tau_E / |M_{<=t}|
for the Permuted task:  tau_E = 1/rho_p,  M_refresh_dot ~ eta*N,  |M| ~ N
=>  c ~ eta / rho_p
```

## Why a naive ablation metric is a quasi-tautology

A tempting but invalid operationalization is `nu^op` = *the fraction of hidden
units whose ablation does not raise the **current-task** loss above a small
threshold*. It is tautological. If drift is fast enough that the network never
learns the current permutation (accuracy near chance — ≈ 0.1 for 10 classes),
there is nothing to lose: ablating any unit cannot raise the loss, so
`nu^op -> 1` **mechanically**. Such a metric measures `(1 − learnedness)`, not
the "fraction of units encoding stale structure", and its growth with `rho_p`
would confirm only the triviality *fast drift → network fails to train*, not
Lemma 2. The accuracy-gate and past-control below exist precisely to avoid this.

## The honest metric

The metric rests on two controls.

1. **Accuracy-gate / dwell.** Nostalgia is memory that is *retained* and was
   *once predictive* — meaningful only when the current task is actually learned.
   The dwell time `dwell = 1/rho_p` is chosen `>>` the steps-to-convergence on a
   single permutation (with `eta = 0.1`, convergence takes ~1000 steps), so most
   scan points are learned (accuracy `>>` chance). Points with
   `accuracy < ACC_GATE = 0.80` are flagged **undertrained** and **excluded**
   from any conclusion about `nu`.

2. **True nostalgia via past-control.** A unit is **nostalgic** iff its ablation
   (a) does **not** raise the current-permutation loss above threshold (useless
   now) **and** (b) **does** raise the loss on at least one **recent past**
   permutation (it encoded structure that was useful then). We keep test batches
   of the `n_past = 3` most recent past permutations.
   - `nu^op_true = mean[ (dL_current ≤ thr) & (max_k dL_past_k > thr) ]`;
   - `dead_frac = mean[ (dL_current ≤ thr) & (max_k dL_past_k ≤ thr) ]`
     (dead/redundant — useful nowhere);
   - `useful_frac = mean[ dL_current > thr ]`.

## Setup

- **Dataset:** `sklearn.datasets.load_digits()` — real handwritten 8×8 digits
  (64 pixels, 10 classes, 1797 samples), shipped *with* scikit-learn, **no
  network download**. An honest mini-analogue of Permuted MNIST. Synthetic
  teacher fallback if sklearn is unavailable.
- **Model:** pure-numpy MLP `64 -> H -> 10` (tanh hidden layer, softmax +
  cross-entropy), online SGD, learning rate `eta`. `H = 32` (large enough that a
  trained network has a non-degenerate stale-unit signal; at `H = 16` a learned
  net has essentially no nostalgic units across the whole learned band).
- **Drift:** a global coordinate permutation `pi_t` shuffles the pixels; labels
  unchanged. With probability `rho_p` per step the permutation is replaced by a
  fresh random one; the displaced permutation enters a queue of recent past
  permutations.

## Parameters

| parameter | value | meaning |
|---|---|---|
| `eta` | `0.1` | SGD learning rate (numerator of `c`); large enough that convergence fits inside the dwell |
| `H` | `32` | hidden units (ablation unit) |
| `N_STEPS` | `120 000` | SGD steps per run |
| `abl_threshold` | `0.05` | `dL` threshold (useless/useful) |
| `n_past` | `3` | recent past permutations kept for the past-control |
| `ACC_GATE` | `0.80` | below it a point is *undertrained*, excluded from `nu` conclusions |
| `N_runs` | `4` | averaging over seeds |
| `rho_p` grid | `{3e-3, 1e-3, 5e-4, 3e-4, 1e-4}` | `dwell = 1/rho_p` from 333 to 10 000 |
| `SEED` | `20260524` | fixed seed (byte-reproducible output) |

## Results

| `rho_p` | `c = eta/rho_p` | `dwell` | accuracy | learned? | `nu^op_true` | dead | useful |
|--------:|----------------:|--------:|---------:|:--------:|-------------:|-----:|-------:|
| 3e-3 | 33.3 | 333 | 0.777 | **no** | 0.007 | 0.000 | 0.993 |
| 1e-3 | 100.0 | 1 000 | 0.880 | yes | 0.025 | 0.003 | 0.972 |
| 5e-4 | 200.0 | 2 000 | 0.912 | yes | 0.035 | 0.006 | 0.959 |
| 3e-4 | 333.3 | 3 333 | 0.955 | yes | 0.141 | 0.026 | 0.833 |
| 1e-4 | 1000.0 | 10 000 | 0.989 | yes | 0.402 | 0.088 | 0.510 |

Over the 4 **learned** points (acc ≥ 0.80):

- Spearman(`rho_p`, `nu^op_true`) = **−1.000** — true nostalgia **decreases**
  with the switch rate (it *grows as drift slows*), the **opposite** of the
  naive Lemma-2 reading.
- Spearman(`c`, `nu^op_true`) = +1.000.
- Spearman(accuracy, `nu^op_true`) = +1.000 — the *trivial* anti-correlation
  with learnedness (characteristic of the naive metric) is absent here; the
  observed correlation has the opposite sign and a different cause (see below).

## Honest framing

**The surrogate does not confirm the naive prediction "faster drift → more
nostalgia". With learnedness controlled, that effect disappears and reverses.**
The monotone `nu ↑ rho_p` that a naive ablation metric would show is the
`(1 − learnedness)` artifact: it is driven by the network failing to learn under
fast drift, not by retained stale structure.

What the surrogate *does* show is a real, well-defined mechanism, but a
**different** one:

- At **slow drift / long dwell** (`rho_p = 1e-4`, dwell = 10 000) the network
  over-specialises to each permutation. After a switch, many units that were
  sharply tuned to the *recent past* permutation are momentarily useless on the
  new one yet still detectably useful on the past — exactly the stale-retained
  signature `nu^op_true` is designed to catch (`nu^op_true ≈ 0.40`).
- At **fast (but still learned) drift** (`rho_p = 1e-3`, dwell = 1 000) the
  network only ever consolidates broadly useful units; little permutation-
  specific structure accumulates, so `nu^op_true ≈ 0.03`.

So the true direction here is *nostalgia grows with the **time spent
over-fitting a regime**, i.e. with dwell `1/rho_p`* — not with the switch rate.
This is informative for § 8.4: the surrogate isolates the genuine stale-retained
fraction (separated from never-useful/dead units) and shows that, in this
small-MLP online-SGD setting, the dominant driver of true nostalgia is
regime-specific over-specialisation under long dwell, **not** memory lagging fast
drift. The naive "monotone-with-`rho_p`" claim does not survive a learnedness
control and should not be reported as confirmation of Lemma 2.

This is an **order-of-magnitude mechanism surrogate**, not a precision
calibration; the small `load_digits` (8×8) MLP is a deliberately minimal probe.

## How to run

```
cd simulations/permuted_mnist
python main.py
```

A rerun yields **byte-identical** numbers (seed `20260524`, runs decorrelated by
`SEED + run`; no wall-clock in the summaries).

## Files

- `model.py` — dataset, numpy MLP, counterfactual ablation, `nu^op_true` with
  past-control, drift online simulation for one `rho_p`.
- `main.py` — `rho_p` scan, accuracy-gate, seed averaging, gated correlation
  check, summary writing, optional figure.
- `results_summary.txt` / `.json` — human- and machine-readable tables.
- `run.log` — copy of the run's stdout.
- `fig_nu_vs_c.png` — `nu^op_true` and accuracy vs `c = eta/rho_p`, with
  undertrained points marked and the accuracy gate drawn.
