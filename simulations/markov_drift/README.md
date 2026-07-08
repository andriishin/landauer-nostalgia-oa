# Markov-drift simulations (§ 6.1 of *The Non-Stationary Predictive Efficiency*)

Numerical illustration of the three regimes of the predictive information
$I_{\mathrm{pred}}(t)$ and nostalgia $\nu^{\mathrm{op}}(t)$ for the PSP (piecewise-stationary with Poisson
switches) surrogate of the OU-drift class (§ 5.2).

## Model in one paragraph

- Finite-state Markov chain on $\Omega=\{1,\dots,k\}$ with $k=8$.
- Transition matrix $P(t)$ is piecewise-constant on intervals
  $[T_i,T_{i+1})$; $\{T_i\}$ is a Poisson process of intensity $\lambda$;
  at each $T_i$ a fresh row-stochastic matrix is sampled from
  Dirichlet$(\alpha)$ on each row (we use $\alpha=0.3$, see below).
- Observer maintains $\hat P(t)$ — sliding-window Laplace-smoothed estimator
  of the transition matrix over the last $\tau_w$ steps; for the "no reset"
  regime $\tau_w=\infty$.
- Predictive information at horizon $\tau$:
  $I_{\mathrm{pred}}(t,\tau) = I_{\mathrm{opt}}(P(t),\tau)
   - \mathbb{E}_{X_t\sim\pi(t)} D_{\mathrm{KL}}(P(t)^\tau(\cdot|X_t)\,\|\,\hat P(t)^\tau(\cdot|X_t))$,
  clipped at zero.
- **Canonical nostalgia** — prediction shortfall
  $\nu^{\mathrm{op}}(t)=1-I_{\mathrm{pred}}(t,\tau)/I_{\mathrm{opt}}(P(t),\tau)$ (leading measure; discriminates regimes).
- **Efficiency**: $\eta = I_{\mathrm{pred}}/I_{\mathrm{mem}}$ with the *information* denominator
  $I_{\mathrm{mem}}(t)=(K/2)\ln(\max(N_{\mathrm{obs}},e))$, $N_{\mathrm{obs}}=t$ — the Bayesian complexity of a
  $K$-parameter model over $N_{\mathrm{obs}}$ observations.
- Structural ballast $\nu^{\mathrm{Still}}=1-\eta$ is secondary (structurally $\approx 1$).

### Why $\alpha=0.3$ and $\tau=1$?

The text of § 6.1 nominates $\alpha=1$ (uniform on the simplex). At $\alpha=1$
Dirichlet rows are very flat, so $\tau_{\mathrm{relax}}\sim 1.5$ steps and
$P^{\tau\ge 2}$ converges to the stationary distribution within numerical
noise. To keep the comparison non-trivial we draw from Dirichlet$(\alpha=0.3)$,
which gives $\tau_{\mathrm{relax}}\sim 2{-}3$ steps (still adiabatic at
$\lambda=10^{-3}$, $\lambda\tau_{\mathrm{relax}}\sim 3\times 10^{-3}\ll 1$),
and report $I_{\mathrm{pred}}$ at one-step horizon $\tau=1$. This is a
methodological choice for the surrogate, not a change to the analytical
class § 5.2.

## How to run

```bash
pip install -r requirements.txt
python main.py
```

The script runs ~250 individual simulations of length $10^4$–$10^5$ steps
each. End-to-end runtime on a modern laptop is roughly 5–15 minutes.

## Outputs

- `../../paper/figs/fig1_Ipred_three_regimes.png` (+ `Fig1.pdf`) — $I_{\mathrm{pred}}(t)$ for the three
  regimes (stationary; drift, no reset — collapse; drift with reset at three $\tau_w$).
- `../../paper/figs/fig2_Ipred_avg_vs_tau_w.png` (+ `Fig2.pdf`) — late-time $\langle I_{\mathrm{pred}}\rangle$
  vs the memory window $\tau_w$ (dome with an interior optimum) at $\lambda=10^{-3}$.
- `../../paper/figs/fig3_nu_evolution.png` (+ `Fig3.pdf`) — $\nu^{\mathrm{op}}(t)$ for the drift-no-reset
  and drift-with-reset regimes.
- `results_summary.txt` and `results_summary.json` — § 6.1 numbers: $I_{\mathrm{pred}}$, $I_{\mathrm{mem}}$,
  $\eta$, $\nu^{\mathrm{op}}$, $\nu^{\mathrm{Still}}$ per regime; $I_{\mathrm{pred}}$ collapse by $\sim 940\times$;
  optimal $\tau_w^\star = 200$ ($\tau_E = 1000$); bound checks.

## Files

- `model.py` — random transition matrices, sliding-window estimator,
  predictive-information functionals, single-run simulator.
- `main.py` — experiment driver, plotting, quantitative analysis.

## Reproducibility

Seeded with `SEED=20260524` and per-condition offsets (see `main.py`).
Re-running on the same NumPy version reproduces the figures and the
numbers in `results_summary.json` byte-identically. A Zenodo deposit
will be created at journal submission.
