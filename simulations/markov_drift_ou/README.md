# OU-drift simulation (§ 6.2 of *The Non-Stationary Predictive Efficiency*)

Numerical illustration of the collapse regime (Lemma 2, § 5.1) on the canonical OU
class of § 5.2; complement to the PSP surrogate in `../markov_drift/`, which only
illustrates the regime qualitatively. At $\sigma = 0.1$ ($\sigma^2/\lambda = 10$) the
OU-concentration condition is **not** satisfied, so conjecture S8.1 / (S8.2)
Supplementary on the adiabatic $(K/2)\ln(\lambda t)$ asymptotic is **not tested here** —
that is done in `../markov_drift_ou_iinf/` and `../markov_drift_ou_iinf_adiab/`.

## Model in one paragraph

- Finite-state Markov chain on $\Omega=\{1,\dots,k\}$ with $k=8$.
- Off-diagonal logits $\theta_{ij}(t)$ ($i\ne j$) of the transition matrix
  evolve by independent Ornstein–Uhlenbeck SDEs
  $d\theta_{ij} = -\lambda(\theta_{ij}-\theta_{ij}^*)\,dt + \sigma\,dW_{ij}$,
  with $\theta_{ij}^*=0$ (uniform stationary environment),
  $\theta_{ii}=0$ as anchor for softmax row-normalisation.
- $P(t) = \mathrm{softmax}_{\text{row}}(\theta(t))$ (over all $k$ columns including
  the anchored diagonal); $K = k(k-1) = 56$ free parameters.
- Euler–Maruyama integration with $\Delta t = 1$.
- Observer maintains $\hat\theta(t)$ via online maximum-likelihood ascent on the
  per-step log-likelihood with decaying step $\eta_t = c/\sqrt{t}$; this is the
  diffusion-prior-Bayes approximation for the OU posterior referenced in § 5.2.
- Predictive information at horizon $\tau=1$:
  $I_{\mathrm{pred}}(t,\tau) = I_{\mathrm{opt}}(P(t),\tau)
   - \mathbb{E}_{X_t\sim\pi(t)} D_{\mathrm{KL}}(P(t)^\tau \,\|\, \hat P(t)^\tau)$,
  clipped at zero.
- **Efficiency**: $\eta = I_{\mathrm{pred}}/I_{\mathrm{mem}}$ with the *information* denominator
  $I_{\mathrm{mem}}(t) = (K/2)\ln(\max(N_{\mathrm{obs}},e))$, $N_{\mathrm{obs}}=t$ — the MDL estimator of growing
  memory (Bayesian complexity of a $K$-parameter model).
- **Canonical nostalgia** — prediction shortfall $\nu^{\mathrm{op}}(t) = 1 - I_{\mathrm{pred}}/I_{\mathrm{opt}}$
  (leading measure; discriminates regimes); structural ballast $\nu^{\mathrm{Still}} = 1 - \eta$ is secondary (structurally $\approx 1$).

## Why these parameter values

- $\lambda = 10^{-3}$, $\sigma = 0.1$: OU stationary variance
  $\sigma^2/(2\lambda) = 5$; this is a *strong fluctuation* regime
  ($\sigma^2/\lambda = 10$), **not** an OU-concentration / adiabatic regime —
  see supplementary § S4.2. We pick it because it is the most charitable
  working point: it keeps the chain genuinely informative
  ($I_{\mathrm{opt}} \sim 0.7$ nats/step), whereas the adiabatic
  $\sigma \sim 3\cdot 10^{-3}$ collapses $I_{\mathrm{opt}}$ into numerical
  noise. The transient-window upper bound is $1/\sigma^2 \approx 10^2$, so
  $T = 2\cdot 10^4$ runs well beyond the strict transient bound; posterior
  concentration is nonetheless the leading effect at this point.
- $T = 2\cdot 10^4$, $50$ independent realisations.
- Fit window $[1000, 15000]$: $t \gg \tau_E = 10^3$ (eliminates start-up).
- Initial logits $\theta(0) \sim \mathcal{N}(0, 0.1^2)$.

## How to run

```bash
pip install -r requirements.txt
python main.py
```

Runtime is roughly 5–15 minutes on a modern laptop (dominated by the per-step
softmax + 8x8 stochastic sampling across $50 \times 2\cdot 10^4 = 10^6$
chain steps).

## Outputs

- `../../paper/figs/fig4_ou_nuop_vs_t.png` (+ `Fig4.pdf`) — **carrier figure**: $\nu^{\mathrm{op}}(t)$ of both
  learners under drift; Robbins–Monro freezes ($\nu^{\mathrm{op}}$ rises toward 1), the constant step holds at a floor.
- `../../paper/figs/fig5_ou_optimal_step.png` (+ `Fig5.pdf`) — **design decides**: late-time constant-step
  $\nu^{\mathrm{op}}$ vs step magnitude at $\sigma=0.1$; interior optimum $s^\star \approx 0.6$.
- `results_summary.txt` and `results_summary.json` — $I_{\mathrm{opt}}$, $I_{\mathrm{mem}}$, late-time
  $I_{\mathrm{pred}}$, $\eta$, $\nu^{\mathrm{op}}$ for both learners; sweeps over $\sigma$ and step magnitude; bound checks.
- `run.log` — full progress log of the simulation run.

## Result

Robbins–Monro **freezes** under drift: $\nu^{\mathrm{op}}$ rises $0.567 \to 0.892$ ($t = 2200 \to 20000$),
$\langle\nu^{\mathrm{op}}\rangle = 0.820$, $\eta \approx 5\cdot10^{-4}$. The constant step $s=0.3$ **tracks**:
$\nu^{\mathrm{op}}$ holds at $\approx 0.33$, $\langle\nu^{\mathrm{op}}\rangle = 0.316$. The step-magnitude sweep
gives an **interior optimum** $s^\star \approx 0.6$ ($\langle\nu^{\mathrm{op}}\rangle = 0.290$), comfortably
beating the frozen Robbins–Monro (0.820) — "forget at the pace of drift". Bounds: $\eta \in [0, 0.008]$,
$\nu^{\mathrm{op}} \in [0, 1]$, zero $I_{\mathrm{pred}} > I_{\mathrm{opt}}$ violations over 50000 points (see
`results_summary.txt`). The S8.1 adiabatic asymptotic is not tested at this point ($\sigma^2/\lambda=10$) —
it is checked in `../markov_drift_ou_iinf_adiab/`.

## Files

- `model.py` — OU step on logits, softmax row-normalisation, online ML logit
  estimator, predictive-information functionals, single-run simulator.
- `main.py` — experiment driver, plotting, quantitative analysis.

## Reproducibility

Seeded with `SEED=20260524` and per-run offsets (see `main.py`). Re-running
on the same NumPy version reproduces the figures and the numbers in
`results_summary.json` byte-identically. A Zenodo deposit will be created
at journal submission.

## Relationship to `../markov_drift/`

The PSP simulation in `../markov_drift/` is the *qualitative* numerical
illustration of § 6.1; it does not produce the BNT-like
$(K/2)\ln(\lambda t)$ term because Poisson regime switches without a
recurrence-dictionary destroy posterior concentration at each switch (see
§ 5.2 caveat and § 6.1 quantitative discussion). The OU simulation here
is the *strict* structural counterpart: drift is continuous, the parameter is the
same $\theta(t)$ throughout, and the ideal Bayesian learner sees a
recurrent latent — i.e. the structural preconditions of S8.1 hold here. At
$\sigma = 0.1$, however, the system is outside the OU-concentration regime, so the
BNT logarithm $(K/2)\ln(\lambda t)$ does not appear in this parametrisation (see the
"Result" section); the test of the S8.1 asymptotic itself is in
`../markov_drift_ou_iinf_adiab/`.
