# iter-006 Step 1: OU + I_infinity (BNT redundancy) — critical experiment

Companion to `../markov_drift_ou/` (iter-005). That simulation showed that the
asymptotic of conjecture S8.1 / (S8.2) Supplementary *fails* with the one-step MI definition of
`I_pred`: the empirical slope was -0.089 vs theoretical 28 (ratio -0.003, both
sign and magnitude wrong).

This directory tests **the structural fix**: redefine `I_pred` as a BNT-style
*cumulative* predictive information (excess entropy / Bayesian redundancy)
rather than a one-step MI, and re-test the slope.

## Two metrics

### Metric A — Growing-window MI (control)

`I^A_pred(t, T_w) = I(X_{t-T_w+1..t} ; X_{t+1..t+T_w})`

For a Markov chain this *analytically* equals `T_w * I_1 + O(1)` where `I_1` is
the one-step MI under the current `P(t)`. We compute it in closed form (the
joint distribution over `8^T_w` symbols is intractable for `T_w > 3`). This is
a **control** that makes the one-step ceiling explicit: `I^A_pred ≤ T_w ln k`,
so the per-symbol information is still capped at `ln k`.

### Metric B — BNT redundancy (main test)

Cumulative excess code-length

```
L_excess(t) := sum_{s=1..t} [ log P_s(X_s | X_{s-1}) − log hat_P_s(X_s | X_{s-1}) ]
```

This is the sample-by-sample analogue of cumulative `KL(P_s || hat_P_s)`. The
classical Bayes-redundancy theorem (BNT 2001, Class II / regular K-dim model)
predicts in a stationary regime

```
E[L_excess(t)]  ~  (K/2) ln t  +  O(1).
```

In the slowly-drifting OU regime `t << 1/sigma^2` the same scaling carries
over with `t → lambda*t`, giving slope **K/2 = 28** in `L_excess(t)` vs
`ln(lambda t)`.

We also report the *operational* form

```
G_obs(t)   := sum_{s=1..t} [ log hat_P_s(X_s | X_{s-1}) + log k ]
eta^inf(t) := G_obs(t) / (R t)
```

so the operational slope test is `eta^inf(t)*t` vs `ln(lambda t)`, theoretical
slope **K/(2R) = 28**. At finite `sigma` the OU drift adds a *linear-in-t* KL
channel on top of the BNT log-term, and a term linear in `t` cannot be absorbed
into the intercept of a fit against `ln(lambda t)` — it contaminates a naive
log-only slope. The BNT coefficient is therefore read off a **joint fit**
`L_excess(t) = A*t + B*ln(lambda t) + C`, with `B` the quantity compared to `K/2`.

## Parameters

| symbol | value | rationale |
|---|---|---|
| k | 8 | same as paper |
| K = k(k-1) | 56 | OU free logits |
| lambda | 1e-3 | OU mean-reversion rate (paper §5.2) |
| sigma | **0.01** | **smaller than iter-005's 0.1**: keeps T<<1/sigma^2 over whole sim |
| T | 1e4 | = 1/sigma^2, exactly the transient bound |
| n_runs | 30 | small per spec (sufficient for slope of order O(10)) |
| theta0_std | 0.1 | initial logit spread |
| measure_every | 100 | subsampling |
| R | 1 | budget rate (unitless) |
| learn_c | 3.0 | Robbins-Monro initial step |
| fit window | [1e3, 1e4] | t >> tau_E and t <= 1/sigma^2 |

## How to run

```bash
pip install -r requirements.txt
python main.py
```

Runtime: ~30-90 s on a laptop (1e4 steps × 30 runs).

## Outputs

- `fig_iinf_growing.png` — Metric A: `I^A(t=T_final, T_w)` vs `T_w`
- `fig_iinf_etat_vs_lnlambdat.png` — Metric B operational fit: `eta^inf*t` vs `ln(lambda t)`
- `fig_iinf_Lexcess.png` — Metric B BNT-strict fit: `L_excess(t)` vs `ln(lambda t)`
- `results_summary.{txt,json}` — fit slopes, R^2, ratios to theory
- `run.log` — full log

## What the script checks

The test is the BNT log-coefficient `B` from the joint fit
`L_excess(t) = A*t + B*ln(lambda t) + C`, compared against `K/2 = 28`. The naive
log-only slope of `L_excess(t)` against `ln(lambda t)` and the operational
`eta^inf(t)*t` slope are also reported, but at finite `sigma` both are dominated
by the linear drift channel and are *not* directly comparable to `K/2`. The ratio
`B/(K/2)` measures how closely the cumulative excess-code-length follows the BNT
Class II log-scaling once the linear channel is separated.

At `sigma = 0.01` (`results_summary.txt`) the joint-fit log-coefficient is
`B = 9.5` (`B/(K/2) = 0.34`) — well below `K/2`. The naive log-only slope is
instead `87.8` (ratio `3.1`, *above* `K/2`) and the operational slope is `-8.7`;
both are dominated by the fitted linear term `A ~ 1.85e-2`, which is exactly why
the joint fit is needed. The diagnosis: the OU drift adds a *structural*
linear-in-t KL loss that masks the BNT log-term at finite `sigma`. The adiabatic
scan in `../markov_drift_ou_iinf_adiab/` separates the two channels and shows
`B/(K/2)` rising toward `1` as `sigma^2/lambda -> 0`.
