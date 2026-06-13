# OU-drift simulation (§ 6.2 of *The Non-Stationary Landauer Efficiency*)

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
- $\nu(t)=1-I_{\mathrm{pred}}/I_{\mathrm{opt}}$.
- Two budgets reported:
  - $N_{\max}^{Rt}(t) = R\cdot t$ (matches theory derivation of (S8.2));
  - $N_{\max}^{\text{full}}(t) = \sum_s [R + \nu(s)\cdot |M|/\tau_{1/2}]$
    with $|M|=K=56$, $\tau_{1/2}=100$, $R=1$.
- $\eta_L(t) = I_{\mathrm{pred}}/N_{\max}(t)$.

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

- `../../paper/figs/fig4_ou_eta_vs_t.png` — $\eta_L(t)$ for OU drift averaged
  across realisations, with $\pm 1\sigma$ band; both $N_{\max}$ variants.
- `../../paper/figs/fig5_ou_etat_vs_lnlambdat.png` — $\eta_L\cdot t$ vs
  $\ln(\lambda t)$ with the empirical linear fit and the theory line
  $K/(2R) = 28$ for visual slope comparison.
- `results_summary.txt` and `results_summary.json` — fit slope, $R^2$,
  $K_{\text{eff}}$, late-time $\nu$, comparison to theoretical $K/(2R)$.
- `run.log` — full progress log of the simulation run.

## Result

The fitted slope of $\eta_L\cdot t$ vs $\ln(\lambda t)$ with $N_{\max}=R\,t$ is
**$-0.089$** against the theoretical $K/(2R)=28$ (ratio $-0.003$; both sign and
magnitude wrong), $R^2\approx 0.89$; $\nu\approx 0.82$ at late times — see
`results_summary.txt`. This is **expected**: the point $\sigma=0.1$ lies outside the
OU-concentration regime in which S8.1 is stated, so the absence of the BNT logarithm
here does **not** refute the conjecture — it only confirms that, away from the
adiabatic limit, $\eta_L(t)$ stays in the collapse regime (Lemma 2). The test of the
S8.1 asymptotic itself is carried out in `../markov_drift_ou_iinf_adiab/`, where the
slope approaches $K/2$ as $\sigma^2/\lambda \to 0$.

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
