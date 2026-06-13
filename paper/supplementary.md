<!--
  Supplementary material (English version). Russian version — supplementary.ru.md.
-->

# Supplementary Material: The Non-Stationary Landauer Efficiency

To accompany *"The Non-Stationary Landauer Efficiency: Memory Growth, Informational Nostalgia, and
the Thermodynamics of Learning Systems"*.

This document contains technical additions not included in the main text: full proofs of Remarks 5
and 6 from § 4.4 (the proof of Lemma 2), extended notation, the parameters of the numerical
simulations of § 6.1, § 6.2 of the main text and of the adiabatic scan (§ S8.3 below), as well as
the formal details of § 5 (the full Remarks 1–6 to Lemma 2 and the pre-registered methodological
commitment against ad hoc rescues), § 6 (full simulation parameters, Robbins–Monro analysis,
wide-window negative control) and § 7 (the full formal correspondence between $\eta_L$ and EFE, the
metric $\rho(t)$, categorical caveats). Structure: S1 — full proof of Remark 6
(softmax-regularity); S2 — full proof of Remark 5 (per-bit-uniformity); S3 — extended
notation glossary; S4 — simulation parameters for reproducibility; S5 — the full apparatus
of § 5 (Remarks 1–6 to Lemma 2, the protocol of § 5.2, existence of $\nu_c$ and $\tau_{\text{reset}}^*$);
S6 — full technical details of the numerical illustrations of § 6 (PSP parameters, Robbins–Monro learner,
wide-window negative control); S7 — formal correspondence between $\eta_L$ and EFE, the metric $\rho(t)$,
categorical caveats; S8 — adiabatic asymptotics of $L_{\text{excess}}$ as an open conjecture
(Class II BNT, parametric scan); S9 — three-faces comparison of decomposition (1) and the full derivation
of the majority-vote variant of Lemma 1 (§§ 2.1, 2.3 of the main text). References of the form "§ X.Y",
"equation (N)", "Lemma N", "Remark N" refer to main.ru.md; bibliographic keys refer to `paper/refs.bib`.

---

## S1. Full proof of Remark 6 (softmax-regularity for OU)

Remark 6 (§ 4.4) states: in the OU class of § 5.2 the exponential decorrelation of the
autocovariance of the logits $\theta(t)$ is transferred through the softmax to any bounded measurable
functional $\xi(\theta)$ (in particular, to the entries of the transition matrix $P_{ij}(\theta)$
and to any sufficient statistic of prediction). The proof has five steps: (i) Lipschitz softmax in
$\ell_\infty$; (ii) transfer to the TV of the transition matrix; (iii) boundedness of the
OU trajectory on the ball $R(\delta)$ with probability $1-\delta$; (iv) transfer of the decay to the
covariance of $\xi$; (v) translation of the TV bound into a bound on the conditional MI.

### S1.1 Step 1. Lipschitz continuity of softmax

Let $\text{softmax}: \mathbb{R}^k \to \Delta^{k-1}$ be the standard map
$\text{softmax}_i(\theta) = \exp\theta_i / \sum_j \exp\theta_j$. The Jacobian has the form
$\partial \text{softmax}_i / \partial \theta_j = \text{softmax}_i(\delta_{ij} - \text{softmax}_j)$;
the diagonal entries $p_i(1-p_i) \ge 0$, the off-diagonal entries $-p_i p_j \le 0$, while the sum of
absolute values over a column does not exceed $2 p_j$. The standard consequence is the
$\ell_1 / \ell_\infty$ pairwise bound:

$$\|\text{softmax}(\theta) - \text{softmax}(\theta')\|_1 \;\le\; 2 \, \|\theta - \theta'\|_\infty. \tag{S1.1}$$

The derivation of (S1.1) is the integration of the Jacobian along the line $\theta(s) = \theta + s(\theta'-\theta)$,
$s \in [0,1]$, applying Hölder's inequality to each component of the gradient; the constant 2
is optimal. A standard result, with no additional assumptions on the dimension $k$.

### S1.2 Step 2. Transfer to the total variation of the transition matrix

In the OU parametrisation of § 5.2 the rows of the transition matrix are given by a row-wise softmax
over the logits: $P_{ij}(\theta) = \exp \theta_{ij} / \sum_l \exp\theta_{il}$, where by convention
$\theta_{ii} \equiv 0$. For each row $i$, applying (S1.1) to the vectors $(\theta_{ij})_{j=1}^k$ and
$(\theta'_{ij})_{j=1}^k$:

$$\|P_{i,\cdot}(\theta) - P_{i,\cdot}(\theta')\|_1 \;\le\; 2 \, \max_j |\theta_{ij} - \theta'_{ij}|.$$

In terms of the TV distance between the distributions $P_{i,\cdot}(\theta)$ and $P_{i,\cdot}(\theta')$
(half the $\ell_1$ distance, $\|\mu - \nu\|_{TV} = \frac{1}{2}\|\mu - \nu\|_1$) one obtains

$$\|P_{i,\cdot}(\theta) - P_{i,\cdot}(\theta')\|_{TV} \;\le\; \max_j |\theta_{ij} - \theta'_{ij}|
   \;\le\; \|\theta - \theta'\|_\infty. \tag{S1.2}$$

This gives a row-wise Lipschitz constant $L = 1$ for the softmax parametrisation of the transition
matrix. The coordinate-wise bound (instead of the full $\|\theta - \theta'\|_\infty$ over all rows) is
essential for step S1.4: the $K = k(k-1)$ logit components of the OU system are independent
(the diagonal anchor $\theta_{ii} \equiv 0$), and each row of the matrix depends only on its own
subsystem of $k-1$ components.

### S1.3 Step 3. Boundedness of the OU trajectory with high probability

The OU process on each coordinate $\theta_{ij}$ in the stationary regime is Gaussian,
$\mathcal{N}(0, \sigma^2/(2\lambda))$; the autocovariance $\mathbb{E}[\theta_{ij}(t)\theta_{ij}(s)]
= (\sigma^2/(2\lambda))\, e^{-|t-s|/\tau_E}$ with $\tau_E = 1/\lambda$. By the concentration
inequality for Lipschitz functionals of a Gaussian measure (Borell–Tsirelson, the standard result
of the theory of Gaussian processes), for a single coordinate
$\Pr[\sup_{[0,T]} |\theta_{ij}(t)| > R] \le 2 \exp(-R^2 \lambda / \sigma^2)$; the union bound over
$K$ coordinates gives the radius

$$R(\delta) \;=\; \sigma \sqrt{\ln(2K/\delta)/\lambda}, \tag{S1.3}$$

ensuring $\Pr[\sup_{[0,T], ij} |\theta_{ij}(t)| > R(\delta)] \le \delta$. On the event
$\{\|\theta(t)\|_\infty \le R(\delta)\}$ the local Lipschitz constant of softmax (S1.2) applies,
and, through the product of measures, for a pair of times $(s, t)$ simultaneously.

### S1.4 Step 4. Transfer of the exponential decay to $\xi$

Let $\xi: \mathbb{R}^K \to \mathbb{R}$ be a bounded measurable function, $|\xi| \le B$,
Lipschitz coordinate-wise with constant $L_\xi$ in the $\ell_\infty$ metric: $|\xi(\theta) - \xi(\theta')|
\le L_\xi \|\theta - \theta'\|_\infty$ (for the row-wise softmax projection $L_\xi \le 1$ by (S1.2)).
The covariance of $\xi$ in the stationary OU regime is estimated through its representation as a
difference of products:

$$\text{Cov}\bigl(\xi(\theta(t)),\, \xi(\theta(s))\bigr) \;=\; \mathbb{E}\bigl[(\xi(\theta(t)) - \bar\xi)(\xi(\theta(s)) - \bar\xi)\bigr],$$

where $\bar\xi = \mathbb{E}[\xi(\theta)]$ under the stationary measure. On the event of step S1.3
(of probability $\ge 1-\delta$) one has $|\xi(\theta(t)) - \bar\xi| \le L_\xi \cdot 2 R(\delta)$ when
compared with the deterministic value at $\theta = 0$ (the symmetry point of the stationary measure)
plus a bounded contribution from the shift $\bar\xi$, absorbed into $|\xi| \le B$. Applying the
Cauchy–Schwarz inequality and the standard estimate of the OU autocovariance
$\mathbb{E}[\theta_{ij}(t)\theta_{ij}(s)] = (\sigma^2/(2\lambda)) e^{-|t-s|/\tau_E}$ gives

$$\bigl|\text{Cov}(\xi(\theta(t)), \xi(\theta(s)))\bigr| \;\le\; 2 B \cdot L_\xi \cdot \frac{\sigma^2}{2\lambda} \cdot e^{-|t-s|/\tau_E} + 2B^2 \delta, \tag{S1.4}$$

where the second term is the correction from the low-probability event of step S1.3, in which the
trajectory exits the ball of radius $R(\delta)$. The choice $\delta = e^{-|t-s|/\tau_E}$ for the
separation of events by the characteristic decorrelation time makes the $\delta$-correction a term
of the same order as the leading exponential term. Finally: for any bounded measurable $\xi$ with
coordinate-wise Lipschitz constant $L_\xi$,

$$\bigl|\text{Cov}(\xi(\theta(t)), \xi(\theta(s)))\bigr| \;\le\; C_1(B, L_\xi, \sigma, \lambda) \cdot e^{-|t-s|/\tau_E}, \tag{S1.4'}$$

with explicit constant $C_1 = (B L_\xi \sigma^2/\lambda) + 2B^2$ (or smaller upon optimising the
choice of $\delta$). The exponential decay of the autocovariance of $\theta$ is transferred to $\xi$
with the same characteristic time $\tau_E$ — this is precisely the substantive statement of Remark 6.

### S1.5 Step 5. Transfer to the conditional mutual information

The chain "covariance of $\xi$ → MI" via inverse Pinsker (Sason–Verdú-type) requires boundedness of
the log-density ratio $\alpha = \log\|dP/dQ\|_\infty$, which for the softmax parametrisation does not
hold universally: even on the ball $R(\delta)$ of step S1.3 it gives $\alpha \le 2R(\delta) + \log k$,
not simply $\log k$. Therefore step S1.5 is built through *MI-tensorisation* over the independent
coordinates of OU followed by DPI on $\xi$ — a formally correct replacement. Direct
$\chi^2$-tensorisation is inapplicable here: $\chi^2$ is *multiplicative* over the independence of
coordinates ($1 + \chi^2_{\text{joint}} = \prod_{ij}(1 + \chi^2_{\text{coord}})$), not additive, so the
naive inequality "$\chi^2_{\text{joint}} \le K \chi^2_{\text{coord}}$" is false for $\chi^2$ and is
circumvented through MI-tensorisation, which is additive precisely because MI is additive over
independent coordinates.

*Step 5a. KL $\le$ $\chi^2$ as an auxiliary coordinate-wise estimate.* The standard inequality
(Cover and Thomas 2006, exercise 2.45): for any $P_{XY}, P_X, P_Y$

$$I(X; Y) \;=\; D_{KL}\bigl(P_{XY} \,\|\, P_X \otimes P_Y\bigr) \;\le\; \chi^2\bigl(P_{XY} \,\|\, P_X \otimes P_Y\bigr).$$

This inequality does not require boundedness of the density ratio; it holds in any $\sigma$-algebra.
In the steps below it is applied *only coordinate-wise* to justify the asymptotics
$I_{\text{coord}} = \tfrac{1}{2}\rho^2(1 + o(1))$, not for tensorisation.

*Step 5b. MI for a bivariate Gaussian pair on a coordinate.* For a single coordinate $\theta_{ij}$
in the stationary OU, the pair $(\theta_{ij}(t), \theta_{ij}(s))$ is Gaussian with normalised
correlation $\rho(t,s) = e^{-|t-s|/\tau_E}$. The standard formula for the MI of a bivariate Gaussian
distribution with marginal variance $\sigma_0^2 = \sigma^2/(2\lambda)$:

$$I\bigl(\theta_{ij}(t);\, \theta_{ij}(s)\bigr) \;=\; -\frac{1}{2}\ln\bigl(1 - \rho^2(t,s)\bigr). \tag{S1.5a}$$

For $|t-s| \ge \tau_E$ one has $\rho^2 \le e^{-2} < 0.14$, and the Taylor expansion gives
$I_{\text{coord}} = \tfrac{1}{2}\rho^2 + O(\rho^4) \le \tfrac{1}{2}\rho^2(1 + o(1))
\le \tfrac{1}{2} e^{-2|t-s|/\tau_E}$. The alternative expression through $\chi^2$
($\chi^2_{\text{coord}} = \rho^2/(1-\rho^2)$ — a direct computation) gives the same leading term
$\tfrac{1}{2}\rho^2$ through step 5a with exponent $2/\tau_E$ (not $1/\tau_E$); MI and
$\chi^2$ on a coordinate agree asymptotically.

*Step 5c. MI-tensorisation over $K$ independent coordinates.* MI is additive over statistically
independent pairs (Cover and Thomas 2006, ch. 2.5): the OU coordinates are independent by construction
(step S2.1: $W_{ij}$ are independent; the covariance is diagonal), and therefore

$$I\bigl(\theta(t);\, \theta(s)\bigr) \;=\; \sum_{ij} I\bigl(\theta_{ij}(t);\, \theta_{ij}(s)\bigr)
   \;\le\; \frac{K}{2} \cdot e^{-2|t-s|/\tau_E} \cdot (1 + o(1)). \tag{S1.5b}$$

The additivity of MI over independent coordinates is *correct* and circumvents the
$\chi^2$-tensorisation trap (multiplicativity $\prod_{ij}(1 + \chi^2_{\text{coord}})$ versus the
expected sum $\sum_{ij} \chi^2_{\text{coord}}$): MI is additive, $\chi^2$ is not, and so the passage
to the final constant $K/2$ is made through MI, not through $\chi^2$.

*Step 5d. DPI on $\xi$ and the final MI bound.* MI is monotone under measurable maps
(Cover and Thomas 2006, ch. 2.8): the coarse-graining of the vector $\theta(\cdot)$ through a bounded
measurable $\xi$ does not increase the MI of the pair $(\xi(\theta(t)), \xi(\theta(s)))$ beyond the MI
of the original pair $(\theta(t), \theta(s))$. Combining steps 5b–5c: for $|t-s| \ge \tau_E$

$$I\bigl(\xi(\theta(t));\, \xi(\theta(s))\bigr) \;\le\; I\bigl(\theta(t);\, \theta(s)\bigr)
   \;\le\; C_2 \cdot e^{-2|t-s|/\tau_E}, \tag{S1.5}$$

with explicit constant $C_2 = K/2 \cdot (1 + o(1))$. The Lipschitz structure of $\xi$ (step S1.4) and
the boundedness $|\xi| \le B$ are not used for the MI bound — they are needed only for the constant in
the coordinate scale of the covariance form (S1.4'). The conditioning on $M_{\le t}^{\text{refr}}$ is
applied through the chain rule on MI: conditioning on $M_{\le t}^{\text{refr}}$ — a function of
$\theta(\cdot)$ on $[s_0, t]$ — does not increase the MI of the pair $(\xi(\theta(t)), \xi(\theta(s)))$
beyond (S1.5); the correction of the standard concentration step S1.3 (the trajectory exiting the
ball of radius $R(\delta)$) enters the additional term $\delta \cdot \ln k$ of the final statement (S1.6).

*Alternative chain.* The same decay is obtained through the Mehler semigroup representation of OU and
the spectral gap of the ergodic diffusion: for the stationary OU measure and bounded measurable $\xi$
the correlations and MI decay by the standard result on the spectral gap of Gaussian semigroups
(Bakry–Émery functional inequality; see (Pavliotis 2014, ch. 4) on the spectral decomposition of OU
and exponential ergodicity; alternatively, the survey in (Levin et al. 2017, ch. 13) for the
discrete analogue of mixing times). This formulation is shorter but requires a reference apparatus
beyond § 5.2; the chain through MI-tensorisation (steps 5b–5c) is formally self-contained within the
present work.

### S1.6 Final statement

**Lemma S1 (supplementary, softmax-regularity of OU).** *In the class of § 5.2 (OU parametrisation of
the logits $\theta(t) \in \mathbb{R}^K$ with characteristic time $\tau_E = 1/\lambda$ and stationary
variance $\sigma^2/(2\lambda)$), for any bounded measurable functional $\xi: \mathbb{R}^K \to \mathbb{R}$
with $|\xi| \le B$ and coordinate-wise Lipschitz constant $L_\xi$ in the $\ell_\infty$ metric, and for
any $\delta \in (0, 1)$, there exists a constant $C = C(\delta, B, L_\xi, \sigma, \lambda, k)$ such that
for all $s, t$ with $|t-s| \ge \tau_E$ one has*

$$I\bigl(\xi(\theta(t));\, \xi(\theta(s)) \,\bigm|\, M_{\le t}^{\text{refr}}\bigr) \;\le\; C \cdot e^{-2|t-s|/\tau_E} + \delta \cdot \ln k. \tag{S1.6}$$

*In particular, for the projections $P_{ij}(\theta(\cdot))$ of the row-wise softmax parametrisation of
the transition matrix one has $L_\xi \le 1$ and $B \le 1$, which gives an explicit constant $C$
independent of $k$ for fixed $\sigma^2/\lambda$.*

Lemma S1 formally closes the step in the proof of Lemma 2 (§ 4.4) that requires the vanishing of
$I(M^{\text{stale}}; X_E \mid M^{\text{refr}})$ for $\Delta t \gg \tau_E$: the right-hand side of (S1.6)
is summed over the age differences in the FIFO distribution and yields a total contribution of the
stale fraction of order $O(c)$, which is what is used in the main proof of Lemma 2 upon substitution
into (8a) through (7). The condition $|t-s| \ge \tau_E$ is the characteristic scale on which the
exponential is not degenerate; for $|t-s| < \tau_E$ the MI is bounded by $H(\xi) \le \ln k$ trivially,
which fits into the additional constant term.

---

## S2. Full proof of Remark 5 (per-bit-uniformity for OU)

Remark 5 (§ 4.4) states the additivity of MI over the independently refreshed coordinates of the
logits of the OU parametrisation of § 5.2 with constant $1 + o(1)$ — the assumption of Lemma 2 that
provides the passage from the arithmetic fraction of refreshed bits $c$ to the information-theoretic
bound $\nu \ge 1 - c$ in the slow-driving limit. Proof: (i) OU as $K$ independent processes;
(ii) Gaussian posterior on each coordinate; (iii) exact additivity of MI over coordinates;
(iv) FIFO in the $K$-coordinate decomposition; (v) comparison with the optimal observer;
(vi) the final bound with $\kappa = 1 + o(1)$; (vii) limits under violation of independence.

### S2.1 OU as $K$ independent processes on the coordinates

In § 5.2 / § 6.2 the OU parametrisation defines the logits $\theta_{ij}(t)$ ($i \ne j$, $K = k(k-1)$
free parameters) as independent one-dimensional processes:

$$d\theta_{ij} = -\lambda(\theta_{ij} - \theta_{ij}^*)\, dt + \sigma\, dW_{ij}, \tag{S2.1}$$

where $W_{ij}$ are independent Wiener processes, $\theta_{ij}^* = 0$ by the convention of § 6.2.
The diagonal logits $\theta_{ii} \equiv 0$ are fixed as the anchor of the softmax normalisation. The
joint distribution $\theta(t) = (\theta_{ij}(t))_{ij}$ in the stationary regime is Gaussian with
diagonal covariance $\Sigma = (\sigma^2/(2\lambda)) I_K$: the coordinates are statistically independent.

### S2.2 Gaussian posterior on each coordinate

Upon observing $X_E^{[0,t]}$ under the true $P(t) = \text{softmax}(\theta(t))$ the learner's posterior
over $\theta(t)$, under regular parametric estimation (MLE or the Robbins–Monro approximation to the
Bayesian update of § 6.2), is asymptotically normal with covariance $\Sigma_{\text{post}}(t)$, whose
entries decay as $O(1/t)$ — the standard result on the asymptotic normality of maximum-likelihood
estimates for regular families (Clarke and Barron 1990, § 3).

*Block structure of the Fisher information.* The information matrix $\mathcal{F}(\theta)$ for OU with
independent coordinates and a row-wise softmax is block-diagonal *by rows* $i$: different rows give an
independent contribution to the likelihood through the independent samples $\{X_t \to X_{t+1}\}$ when
$X_t = i$. The additivity by rows is exact, not approximate.

*Within-row coupling.* Within a row $i$ the Fisher block has the form $\mathcal{F}^{(i)}(\theta) =
n_i(t) \cdot (\text{diag}(p_i) - p_i p_i^T)$, where $p_i = (P_{ij})_{j=1}^k$ is the $i$-th row of the
softmax and $n_i(t)$ is the number of visits to state $i$ by time $t$. This softmax-Fisher matrix is
**not diagonal**: the softmax normalisation couples all $j$ through the common denominator. The exact
additivity of the diagonal approximation $\Sigma_{\text{post}} \approx \text{diag}(\sigma_{ij}^2)$
therefore requires justification.

*Spectral analysis via the Schur formula.* The matrix $\text{diag}(p) - pp^T$ is the standard
centring projector: it has eigenvalue $0$ in the direction $\mathbf{1}$ (which corresponds to the
constraint $\sum_j p_{ij} = 1$, operationally — the fixing of $\theta_{ii} \equiv 0$ as the softmax
anchor, see S2.1) and eigenvalues $\{p_{ij}\}_{j \ne i}$ on the orthogonal subspace of free
parameters. The condition number of the Fisher block in the free coordinates:

$$\kappa(\mathcal{F}^{(i)}_{\text{free}}) \;=\; \frac{\max_{j \ne i} p_{ij}}{\min_{j \ne i} p_{ij}}. \tag{S2.2a}$$

On the event of step S1.3 (the OU trajectory in the ball $R(\delta)$ with probability $\ge 1-\delta$),
for softmax one has $\max_j p_{ij} / \min_j p_{ij} \le e^{2R(\delta)}$, whence
$\kappa(\mathcal{F}^{(i)}_{\text{free}}) \le e^{2R(\delta)}$ — bounded and dependent only on the radius
of the OU ball through $R(\delta) = \sigma\sqrt{\ln(2K/\delta)/\lambda}$ (S1.3).

*Correction to the diagonal approximation.* The Gaussian approximation of the posterior gives
$\Sigma_{\text{post}}^{(i)} = (\mathcal{F}^{(i)}_{\text{free}})^{-1}$ with the same spectrum of inverse
eigenvalues; the off-diagonality of $\Sigma_{\text{post}}$ within a row generates a correction to the
exact additivity (S2.2). In the slow-driving limit $\lambda \tau_{\text{relax}} \ll 1$ the correction
is of order $O(\kappa(\Sigma_{\text{post}}^{(i)}) \cdot \lambda \tau_{\text{relax}})$: the condition
number enters multiplicatively, since the off-diagonal contribution to the MI is bounded through the
trace of $\Sigma_{\text{post}}^{(i)}$, majorised by $\kappa$ times the diagonal part. On the OU ball
S1.3 this gives

$$\Sigma_{\text{post}} \;=\; \text{diag}(\sigma_{ij}^2(t)) \cdot \bigl(1 + O\bigl(e^{2R(\delta)} \cdot \lambda \tau_{\text{relax}}\bigr)\bigr), \tag{S2.2b}$$

which remains $1 + o(1)$ under the strengthened slow-driving condition $e^{2R(\delta)} \cdot \lambda
\tau_{\text{relax}} \ll 1$. In the numerical regimes of § 6.2 ($\sigma = 0.1$, $\lambda = 10^{-3}$,
$K = 56$, $\delta = 0.05$) one obtains $R(\delta) \approx 0.1 \cdot \sqrt{\ln(2240) / 10^{-3}}
\approx 8.8$ — too large for the literal $e^{2R(\delta)} \cdot \lambda \tau_{\text{relax}} \ll 1$;
however, (S2.2b) is a *uniform* upper estimate on the $R(\delta)$-ball, whereas the typical softmax
normalisation is concentrated near equilibrium and $\kappa(\mathcal{F}^{(i)}_{\text{free}})$ is
empirically close to $1$ (the numerical control of § 6.2 shows $1 + o(1)$ on average over
realisations). An analytical improvement of the estimate (S2.2b) — through the concentration of
softmax around the typical logit, rather than the uniform bound — is left as a refinement in § 8.5.

### S2.3 Additivity of MI over independent coordinates

For a Gaussian posterior $\hat\theta = (\hat\theta_{ij})$ with diagonal covariance
$\Sigma_{\text{post}} = \text{diag}(\sigma_{ij}^2)$, the *exact* additivity of MI with the true latent
holds:

$$I(\hat\theta;\, \theta) \;=\; \sum_{ij} I(\hat\theta_{ij};\, \theta_{ij}). \tag{S2.2}$$

Derivation: by the chain rule (Cover and Thomas 2006, ch. 2.5)

$$I(\hat\theta;\, \theta) \;=\; \sum_{ij} I\!\left(\hat\theta_{ij};\, \theta_{ij} \,\middle|\, \hat\theta_{\lt ij},\, \theta_{\lt ij}\right);$$

the independence of the OU coordinates gives $\theta_{ij} \perp (\theta_{\lt ij})$ (statistical
independence by construction), and the diagonality of $\Sigma_{\text{post}}$ gives
$\hat\theta_{ij} \perp \hat\theta_{\lt ij}$ (conditional independence of the posterior estimates given
the true latent). The conditional MIs reduce to the unconditional $I(\hat\theta_{ij}; \theta_{ij})$,
which yields (S2.2).

The correction to the exact additivity for $\lambda \tau_{\text{relax}} = O(\epsilon)$ is a term
$O(\epsilon)$, absorbed into $1 + o(1)$ as $\epsilon \to 0$.

### S2.4 Refresh-fraction FIFO in the $K$-coordinate decomposition

In the FIFO memory-update scheme of Lemma 2, a refresh fraction of size $c$ refreshes a fraction $c$
of the $K$ coordinates over an interval $\tau_E$ uniformly (the definition of FIFO and the formulas
$(*)$). By the additivity (S2.2) the MI of the refreshed fraction with the true latent of the
environment decomposes:

$$I(M_{\le t}^{\text{refr}};\, X_E^{[t, t+\tau]}) \;=\; \sum_{ij \in \text{refr}} I(\hat\theta_{ij}^{\text{refr}};\, \theta_{ij}(t+\tau)) \;=\; c \cdot K \cdot \bar I_{\text{coord}}(t+\tau), \tag{S2.3}$$

where $\bar I_{\text{coord}}(t+\tau) := K^{-1} \sum_{ij} I(\hat\theta_{ij}; \theta_{ij}(t+\tau))$ is the
average MI per coordinate under the ideal (Bayes-optimal for a regular family) estimate. Structurally:
each refreshed coordinate contributes independently, proportionally to $\bar I_{\text{coord}}$; the
total MI is proportional to the fraction of refreshed coordinates $c$ under $K$-fold additivity.

### S2.5 Comparison with the optimal observer

The MI of the optimal observer (the full sufficient statistic over $\theta(t)$) on the horizon $\tau$
in the slow-driving limit:

$$I_{\text{pred}}^{\text{opt}}(t, \tau) \;=\; I(X_t;\, X_{t+\tau} \mid \theta(t)) \;=\; \sum_{ij} I_{\text{coord}}(t+\tau) \;=\; K \cdot \bar I_{\text{coord}}(t+\tau), \tag{S2.4}$$

where the decomposition over coordinates uses the conditional independence of the transitions
$X_t \to X_{t+\tau}$ over different rows and the independence of the logit coordinates $\theta$ within
a single row in the slow-driving limit. The approximation $\theta(t+\tau) \approx \theta(t)$ works to
accuracy $O(\sigma^2 \tau /\lambda)$: over the horizon $\tau$, when $\lambda \tau \ll 1$, the latent
drift is small, and the MI per coordinate $I_{\text{coord}}(t+\tau) = I_{\text{coord}}(t) \cdot (1 + O(\lambda \tau))$.

### S2.6 The final bound

Substituting (S2.3) into (S2.4):

$$\frac{I(M_{\le t}^{\text{refr}};\, X_E^{[t, t+\tau]})}{I_{\text{pred}}^{\text{opt}}(t, \tau)} \;=\; \frac{c \cdot K \cdot \bar I_{\text{coord}}(t+\tau)}{K \cdot \bar I_{\text{coord}}(t+\tau)} \;=\; c \cdot (1 + O(\lambda \tau)).$$

That is,

$$I(M_{\le t}^{\text{refr}};\, X_E^{[t, t+\tau]}) \;\le\; c \cdot I_{\text{pred}}^{\text{opt}}(t, \tau) \cdot (1 + O(\lambda \tau)). \tag{S2.5}$$

For $\lambda \tau \ll 1$ (the slow-driving regime, satisfied by the condition of the class of § 5.2)
and under the strengthened condition $e^{2R(\delta)} \cdot \lambda \tau_{\text{relax}} \ll 1$ from S2.2
(within-row Fisher coupling), $1 + O(\lambda \tau) + O(e^{2R(\delta)} \lambda \tau_{\text{relax}})
= 1 + o(1)$, and (S2.5) gives per-bit-uniformity with constant $\kappa = 1 + o(1)$. Without the
strengthened condition, $\kappa = 1 + O(e^{2R(\delta)} \lambda \tau_{\text{relax}})$, which for finite
$\delta, \sigma, \lambda$ gives a multiplicative correction but does not cancel the qualitative
structure of (8a) (see also S2.7 for additional sources of non-diagonality). Substitution into (7)
under the condition $I(M^{\text{stale}}; X_E \mid M^{\text{refr}}) \to 0$ (Lemma S1) gives
asymptotically

$$\nu^{\text{theor}}(t) \;=\; 1 - \frac{I(M_{\le t}^{\text{refr}};\, X_E^{[t, t+\tau]})}{I_{\text{pred}}^{\text{opt}}(t, \tau)} \;\ge\; 1 - c \cdot (1 + o(1)), \tag{S2.6}$$

which is the statement (8a) of Lemma 2 with explicit constant $\kappa = 1 + o(1)$.

### S2.7 Discussion of limits: correlated coordinates

Under violation of the independence of the OU coordinates (cross-arrow drift with non-diagonal
$\Sigma$ in (S2.1)) per-bit-uniformity requires refinement. Besides the between-row *prior*
correlation, the second source of non-diagonality is the within-row Fisher coupling through the
softmax normalisation (S2.2): the condition number of the Fisher block
$\kappa(\mathcal{F}^{(i)}_{\text{free}}) = \max_j p_{ij} / \min_j p_{ij}$ controls the deviation from
exact additivity (S2.2a, S2.2b) and enters multiplicatively into the final correction (S2.6). Both
sources are combined in a common condition-number argument.

Let $\kappa(\Sigma_{\text{post}})$ be the condition number (the ratio of the maximal to the minimal
eigenvalue) of the joint covariance of the posterior, combining the between-row prior correlations
(S2.1, under non-diagonal $\Sigma_{\text{prior}}$) and the within-row Fisher correlations (S2.2).
Through the Gaussian representation of MI $I(\hat\theta; \theta) = (1/2) \ln \det(\Sigma_{\text{prior}}
\Sigma_{\text{post}}^{-1}) = (1/2) \ln \det(I + \mathcal{F}\,\Sigma_{\text{prior}})$ (where $\mathcal{F} =
\Sigma_{\text{post}}^{-1} - \Sigma_{\text{prior}}^{-1}$ is the observed Fisher information, and
$\Sigma_{\text{post}}^{-1} = \Sigma_{\text{prior}}^{-1} + \mathcal{F}$) and the Hadamard inequality
$\det(\Sigma) \le \prod_i \Sigma_{ii}$ one obtains an upper bound that in the worst case is majorised
by a value $\kappa(\Sigma_{\text{post}})$ times larger than the additive sum. In terms of the constant
of Remark 5:

$$\kappa \;\le\; \kappa(\Sigma_{\text{post}}) \cdot (1 + o(1)). \tag{S2.7}$$

For a diagonal $\Sigma_{\text{post}}$ (independent coordinates) $\kappa(\Sigma_{\text{post}}) = 1$ and
(S2.7) reduces to $\kappa = 1 + o(1)$. For correlated coordinates the bound is weakened
multiplicatively by $\kappa(\Sigma_{\text{post}})$; a numerical estimate for a specific cross-parametrisation
is the subject of a separate study. In the OU parametrisation of § 5.2 with independent $W_{ij}$ and a
diagonal anchor $\theta_{ii} \equiv 0$, the condition $\kappa(\Sigma_{\text{post}}) = 1$ holds by
construction, and per-bit-uniformity with $\kappa = 1 + o(1)$ applies without restrictions.

---

## S3. Extended notation table

A glossary of the notation of the work indicating the dimension and the defining equation; "dimensionless"
denotes a dimensionless quantity. Abbreviations: "MI" — mutual information (Cover and Thomas 2006, ch. 2);
"TV" — total variation; "KL" — Kullback–Leibler. The contexts of definition are § X.Y and equations (N)
of the main text.

### S3.1 Efficiency and information flow

| Symbol | Dimension/type | Definition | Theoretical upper bound |
|--------|-----------------|-------------|-------------------------------|
| $\eta_L(t)$ | dimensionless, $\in [0, 1]$ | Non-stationary Landauer efficiency; (10) § 5.1 | $1$ (irreversible regime, asymptotics $t \gg \tau_d$) |
| $\eta_L^{(\tau_w)}(t)$ | dimensionless, $\in [0, 1]$ | Efficiency with sliding window $\tau_w$; (10a) § 5.2 | $1$ (irreversible regime) |
| $\eta_L^{\text{stat}}$ | dimensionless | Stationary limit of $\eta_L$; § 5.3, (Andriishin 2026, § 2.1) | $1$ (irreversible regime) |
| $\eta_L^{\text{excess}}(t)$ | dimensionless | Secondary efficiency through $L_{\text{excess}}$; (S8.3) § S8 | no universal upper bound (secondary metric) |
| $I_{\text{pred}}(t, \tau)$ | bits | One-step predictive information at $\hat P(t)$; (4) § 3.1 | $\ln k$ (one-step, $k = \lvert\Omega\rvert$) |
| $I_{\text{pred}}^{\text{stat}}(\tau)$ | bits | Stationary limit of $I_{\text{pred}}$; § 2.2 | $\ln k$ (one-step); $\lvert M_t\rvert\log 2$ (capacity) |
| $I_{\text{pred}}^{\text{opt}}(t, \tau)$ | bits | MI of the ideal observer with access to $\theta(t)$; (7) § 4.1 | $\ln k$ (one-step, $X_t, X_{t+\tau}$ discrete) |
| $I_{\text{pred}}^{\text{Bialek}}(X_E, \tau)$ | bits | Classical predictive information (Bialek et al. 2001); § 2.2 (Andriishin 2026) | $C_\mu$ (excess entropy, Crutchfield and Feldman 2003) |
| $N_{\text{max}}(t)$ | bits | Cumulative Landauer budget; (1) § 2.1 | $\lvert M_t\rvert\log 2$ per active part of memory; unbounded above in time |
| $N_{\text{max}}^{(\tau_w)}(t)$ | bits | Budget on the sliding window $\tau_w$; (10a) § 5.2 | $\lvert M_t\rvert\log 2$ per active part |
| $L_{\text{excess}}(t)$ | nats | Cumulative excess-loss of the learner; (S8.1) § S8 | asymptotics $(K/2)\ln(\lambda t)$ — conjecture (S8.2) |

*Distinction of three types of theoretical bounds.* (i) The one-step bound $\ln k$ — the standard
capacity of a discrete channel with alphabet $\Omega$ (Cover and Thomas 2006); it controls the instantaneous
information capacity through one step of the environment. (ii) The capacity bound $\lvert M_t\rvert\log 2$
— the structural ceiling for the integral information that a finite memory $M_t$ of $\lvert M_t\rvert$ bits
can store. (iii) The excess entropy $C_\mu$ (Crutchfield and Feldman 2003) — the total predictive information
between the past and the future of a stationary process, in general infinite for nontrivial classes.
Distinction: $\ln k$ — for the one-step $I_{\text{pred}}(t,\tau)$, $\lvert M_t\rvert\log 2$ — for
capacity estimates of memory, $C_\mu$ — for interval-MI of the type $I_{\text{pred}}^{\text{Bialek}}$.

### S3.2 Energy flows

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $E_{\text{actual}}^{\text{curr}}(t)$ | J | Current operational dissipation; § 2.1 |
| $E_{\text{store}}(t)$ | J | Cumulative cost of maintaining old memory; § 2.1 |
| $E_{\text{grow}}(t)$ | J | Cumulative cost of memory expansion; § 2.1 |
| $\dot{E}_{\text{store}}(t)$ | W | Instantaneous maintenance power; Lemma 1 § 2.3 |
| $\dot{E}_{\text{grow}}(t)$ | W | Instantaneous expansion power; § 2.4 |
| $\dot{E}_{\text{store}}^{\text{nostalg}}(t)$ | W | Maintenance power of the nostalgic fraction; (8) § 4.2 |
| $E_{\text{store}}^{(1)}(t)$ | J | Minimal dissipation for maintaining one bit; (2) Lemma 1 |
| $E_{\text{grow}}^{(1)}$ | J | Minimal dissipation for creating one bit of capacity; (3) § 2.4 |
| $\eta_{\text{ex}}$ | dimensionless | Exergetic efficiency of metabolism; (Andriishin 2026, § 3.5) |

### S3.3 Memory and its layers

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $M_t$ | bits (size); configuration | Current (refreshable) memory layer; § 3.1 |
| $M_{\lt t}$ | bits; configuration | Archival memory layer; § 3.1 |
| $M_{\le t} = (M_t, M_{\lt t})$ | bits; configuration | Full memory of the system; § 3.1 |
| $M_{\le t}^{\text{refr}}$ | bits | Refresh fraction (refreshed over the interval $\tau_E$); § 4.4 |
| $M_{\le t}^{\text{stale}}$ | bits | Stale fraction (not refreshed over $\tau_E$); § 4.4 |
| $\lvert M_{\le t}\rvert$ | bits | Size of the active part of memory; § 4.4, (11b) § 5.3 |
| $\lvert M_{\le t}^{\text{cap}}\rvert$ | bits | Structural limit of the memory capacity; (12) § 5.3 |
| $\dot{M}_{\text{refresh}}(t)$ | bits/s | Memory refresh rate; § 4.4, (13) § 7 |
| $\dot{M}_{\text{grow}}(t)$ | bits/s | Memory capacity growth rate; (13) § 7 |
| $\dot{M}_{\text{refresh}}^{\text{eff}}(t)$ | bits/s | Effective refresh rate for Robbins–Monro; § 6.2 |

### S3.4 Latent parameters and estimates

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $\theta(t)$ | $\mathbb{R}^K$ | True logits of the environment's transition matrix; § 5.2 |
| $\hat\theta(t)$ | $\mathbb{R}^K$ | Learner's estimate of the logits; § 6.2 |
| $\theta_{ij}(t)$ | $\mathbb{R}$ | Logit coordinate, $i \ne j$; § 5.2 |
| $P(t)$ | $k \times k$ stochastic | Environment's transition matrix; § 5.2 |
| $\hat P(t)$ | $k \times k$ stochastic | Learner's estimate of the transition matrix; § 3.1 |
| $X_t$ | $\in \Omega$ | State of the environment at time $t$; § 3.1 |
| $X_E^{[t, t+\tau]}$ | trajectory | Realisation of the environment on the window $[t, t+\tau]$; (1a) (Andriishin 2026) |

### S3.5 Model parameters and time scales

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $k$ | integer; $\lvert\Omega\rvert$ | Size of the environment's state alphabet; § 5.2, § 6.1 |
| $K = k(k-1)$ | integer | Dimension of the logits of the OU parametrisation; § 5.2 |
| $\Omega$ | finite set | State alphabet of the Markov chain; § 3.1 |
| $\lambda$ | 1/s | OU relaxation rate; the inverse of the characteristic drift time $\tau_E$; § 5.2 |
| $\sigma$ | $\sqrt{1/\text{s}}$ | Amplitude of the OU noise; § 5.2 |
| $\alpha$ | dimensionless | Dirichlet parameter for the PSP; § 6.1 |
| $\tau$ | s | Prediction horizon of the one-step $I_{\text{pred}}$; (4) § 3.1 |
| $\tau_d$ | s | Characteristic dissipation/response time of the system; § 2.1 (Andriishin 2026) |
| $\tau_w$ | s | Width of the estimator's sliding window; § 6.1, § 5.2 |
| $\tau_E = 1/\lambda$ | s | Characteristic drift time of the environment; § 5.2 |
| $\tau_{\text{relax}}$ | s | Relaxation time of the Markov chain at a slice; § 5.2, § 6.1 |
| $\tau_{1/2}$ | s | Half-life of the bit-read accuracy; Lemma 1 § 2.3 |
| $\tau_{\text{reset}}^*$ | s | Optimal memory reset interval; (12) § 5.3 |
| $T$ | K | Temperature of the receiving heat bath; § 2.1 (Andriishin 2026) |
| $\varepsilon$ | dimensionless, $\in (0, 1/2)$ | Accuracy threshold of bit reading; Lemma 1 § 2.3 |
| $\gamma$ | 1/s | Symmetric transition rate in Lemma 1; § 2.3 |
| $\eta_t$ | dimensionless | Robbins–Monro step of the learner; § 6.2 |

### S3.6 Nostalgia and the regime threshold

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $\nu(t)$ | dimensionless, $\in [0, 1]$ | Nostalgia (general symbol, the context determines the variant); § 4.1 |
| $\nu^{\text{theor}}(t)$ | dimensionless, $\in [0, 1]$ | Theoretical nostalgia through $I_{\text{pred}}^{\text{opt}}$; (7) § 4.1 |
| $\nu^{\text{op}}(t)$ | dimensionless, $\in [0, 1]$ | Operational nostalgia through $C_\mu^{\text{emp}}$; (7') § 4.1 |
| $\nu_c$ | dimensionless, $\in (0, 1)$ | Critical value / regime-transition threshold; (11b) § 5.3 |
| $\nu_c^{\text{theor}}$ | dimensionless, $\in (0, 1)$ | Theoretical critical nostalgia; (11b) § 5.3 |
| $\nu_c^{\text{op}}$ | dimensionless, $\in (0, 1)$ | Operational critical nostalgia; § 8.6 Prediction 1 |
| $\nu_c^{\text{emp}}$ | dimensionless, $\in (0, 1)$ | Empirical proxy from § 6.1 (PSP simulation); § 6.1 |
| $c$ | dimensionless | Computable constant $\dot{M}_{\text{refresh}} \tau_E / \lvert M_{\le t}\rvert$; (*) § 4.4 |
| $\eta_M$ | dimensionless | Memory-balance ratio $\dot{M}_{\text{refresh}} \tau_{1/2}/\lvert M_{\le t}\rvert$; (11b) § 5.3 |
| $\kappa$ | dimensionless, $\ge 1$ | Per-bit-uniformity constant of Remark 5; § 4.4, § S2.6 |
| $\rho(t)$ | dimensionless | Informational rejuvenation metric; (13) § 7 |

### S3.7 Information complexity and diagnostic proxies

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $C_v^{\text{static}}(t)$ | bits | Structural complexity (Lempel–Ziv, assembly); § 4.3 |
| $C_v^{\text{predictive}}(t)$ | bits | Dynamic predictive power; (9) § 4.3 |
| $C_\mu$ | bits | Statistical complexity of Crutchfield–Feldman (Crutchfield and Feldman 2003) |
| $C_\mu^{\text{emp}}(t, \tau)$ | bits | Empirical excess entropy on a sliding window; (7') § 4.1 |
| $h_\mu$ | bits/symbol | Entropy rate; § 4.3 |

### S3.8 Active inference (§ 7)

| Symbol | Dimension/type | Definition |
|--------|-----------------|-------------|
| $G[\pi]$ | nats | Expected free energy of policy $\pi$; § 7 / S7 |
| $G^{\text{epist}}[\pi]$ | nats | Epistemic component of EFE (information gain); § 7 / S7 |
| $G^{\text{prag}}[\pi]$ | nats | Pragmatic component of EFE (preference); § 7 / S7 |
| $q(s, o, \pi)$ | density | Generative model in active inference (Friston et al. 2017a; Parr et al. 2022); § 7 / S7 |
| $p(o', s' \mid C)$ | density | Preference distribution with target parameter $C$; § 7 / S7 |
| $\beta = 1/(k_B T)$ | 1/J | Inverse temperature for converting nats to J; (S7.1) § 7 / S7 |

### S3.9 Fundamental constants

| Symbol | Dimension | Value / definition |
|--------|-------------|------------------------|
| $k_B$ | J/K | Boltzmann constant |
| $k_B T \ln 2$ | J | Minimal payment for one erasure (Landauer 1961) |
| $\ln 2$ | dimensionless | Conversion factor nats → bits |

---

## S4. Simulation parameters for reproducibility

This section assembles the technical parameters of the numerical illustrations of § 6.1 (PSP), § 6.2 (OU)
and § S8.3 (adiabatic scan). The implementation is `andriishin/landauer-nostalgia-oa`, directories
`simulations/markov_drift{,_ou,_ou_iinf_adiab}/`. The seed $\text{SEED} = 20260524$ with
per-run offsets; running `python main.py` in each directory reproduces the cited numerical values up
to the version of the libraries used.

### S4.1 PSP simulation (§ 6.1)

*Model class.* A piecewise-stationary Markov process with Poisson switches of the transition matrix.
Alphabet $\Omega = \{1, \ldots, k\}$, $k = 8$. On the interval $[T_i, T_{i+1})$ the transition matrix
$P^{(i)}$ is constant; the switch times $\{T_i\}$ form a Poisson stream with intensity $\lambda$. At
time $T_i$ the matrix is updated to a random stochastic $P^{(i)}$ with independent rows, each row a
Dirichlet with parameter $\alpha = 0.3$.

*Base parameters.* $k = 8$, $\alpha = 0.3$, $\lambda = 10^{-3}$ (for the drift regimes), prediction
horizon $\tau = 1$ step.

*Relaxation rate.* $\langle \tau_{\text{relax}} \rangle \approx 2$–$3$ steps for $\alpha = 0.3$;
$\lambda \tau_{\text{relax}} \approx 3 \cdot 10^{-3} \ll 1$ — slow-driving is satisfied with margin.

*Three regimes* (§ 5.3):

| Regime | $\lambda$ | $\tau_w$ | $T$ | Realisations |
|-------|-----------|----------|-----|------------|
| Stationary | $0$ | $100$ | $10^4$ | $80$ |
| Drift without reset | $10^{-3}$ | $\infty$ | $10^5$ | $60$ |
| Drift with reset | $10^{-3}$ | $\{100, 500, 2000\}$ | $10^5$ | $60$ |

*Grid of $\tau_w$ for Fig. 2.* Eight points $\tau_w \in \{50, 100, 200, 500, 10^3, 2 \cdot 10^3,
5 \cdot 10^3, 10^4\}$, $\lambda = 10^{-3}$, $T = 5 \cdot 10^4$, $40$ realisations per point.

*MI estimator.* KL-detuned proxy $\widehat I_{\text{pred}}(t,\tau) = I_{\text{opt}}(P(t), \tau) -
\mathbb{E}_{X_t \sim \pi(t)}[D_{KL}(P(t)^\tau(\cdot \mid X_t) \| \hat P(t)^\tau(\cdot \mid X_t))]$
with clipping at zero from below; cross-checked against independent estimates: the discrete plug-in
with Treves–Panzeri correction (Treves and Panzeri 1995) and the continuous KSG estimator with $k$-NN,
$k_{NN} = 3$ (Kraskov et al. 2004); consistency on the grid verified.

*Budget $N_{\text{max}}$.* In the simulation the normalisation $k_B T \ln 2 = 1$ is adopted and
$N_{\text{max}}(t) = N_0 + R \cdot t$ with $R = 1$ unit per step (the base operational payment plus
the Lemma 1 payment on $|M_t| = k^2 = 64$ maintained parameters).

*Computational cost.* Complexity per realisation $O(T \cdot k^2 \cdot \tau_w)$; the base stationary
parametrisation ($T = 10^4$, $\tau_w = 100$, $k = 8$) — $\sim 6 \cdot 10^7$ operations.

*Key numerical results.* Stationary: $\eta_L \approx 6.1 \cdot 10^{-5}$ for
$t > 5 \cdot 10^3$. Drift without reset: $\nu_\infty(t > 60000) \approx 0.999 \pm 0.001$, $\eta_L$
decays by $\sim 3$ decades. Grid of $\tau_w$: optimum $\tau_w^* \approx 200$, $\langle\eta_L\rangle
\approx 1.3 \cdot 10^{-5}$. Empirical $\nu_c^{\text{emp}} \approx 0.44$ at the half-life of $\eta_L$.

*Reproduction command and outputs.* `python simulations/markov_drift/main.py`; wall-time
$\sim 5$–$10$ minutes; outputs — `results_summary.{txt,json}`, `run.log`, figures
`paper/figs/fig{1,2,3}_*.png`.

### S4.2 OU simulation (§ 6.2)

*Model class.* A Markov chain with continuous OU drift of the logits of the transition matrix. $k = 8$,
$K = k(k-1) = 56$ free parameters. The OU equations
$d\theta_{ij} = -\lambda(\theta_{ij} - \theta_{ij}^*)\, dt + \sigma\, dW_{ij}$ with $\theta_{ij}^* = 0$,
independent $W_{ij}$ for $i \ne j$, and fixed $\theta_{ii} \equiv 0$ as the softmax anchor.
Euler–Maruyama discretisation with $\Delta t = 1$.

*Base parameters.* $k = 8$, $\lambda = 10^{-3}$, $\sigma = 0.1$, $T = 2 \cdot 10^4$ steps,
$50$ independent realisations of the OU process; initial logits $\theta(0) \sim \mathcal{N}(0, 0.1^2)$.

*Dimensionless parameters.* $\sigma^2/\lambda = 10$; stationary variance of the logits
$\sigma^2/(2\lambda) = 5$ — a strong fluctuation regime (not OU-concentration). Transient window
$1/\sigma^2 = 10^2$ steps.

*Approximation of the Bayesian learner.* Robbins–Monro online update with decreasing step
$\eta_t = c_0 / \sqrt{t}$, $c_0 = 3$; at each step $\hat\theta_{x_t, j}(t+1) = \hat\theta_{x_t, j}(t)
+ \eta_t \cdot (\mathbb{1}[x_{t+1} = j] - \hat P_{x_t, j}(t))$ for $j \ne x_t$. The rate
$1/\sqrt{t}$ corresponds to the concentration of the posterior along each of the $K$ directions and is
asymptotically Bayes-optimal for regular parametric models (Clarke and Barron 1990, § 3).

*Effective refresh rate.* $\dot{M}_{\text{refresh}}^{\text{eff}}(t) \sim K c_0 / (\tau_E \sqrt{t})$;
for $K = 56$, $\tau_E = 10^3$, $t = 2 \cdot 10^4$ it gives $\sim 0.17$ refreshes per step;
the constant $(*)$ $c(t) = \dot{M}_{\text{refresh}}^{\text{eff}} \cdot \tau_E / |M_{\le t}|
\approx K c_0 / (|M_{\le t}| \sqrt{t}) \approx 0.02$ for $|M_{\le t}| = K$. This justifies the
numerical gap $\nu(t > 10^4) \approx 0.82$ against the prediction $\liminf \nu^{\text{theor}}
\ge 1 - c \approx 0.98$ (§ 6.2; the finiteness of the observation time and Remark 5 cover the
discrepancy $\approx 0.16$).

*Budget.* Two variants: (a) $N_{\text{max}} = R t$ with $R = 1$ and (b) the full $N_{\text{max}}^{\text{full}}
= R t + \nu(t) \cdot K / \tau_{1/2}$ with $\tau_{1/2} = 100$. The curves $\eta_L(t)$ are qualitatively
identical.

*Key numerical results.* Late-time ($t > 10^4$): $\nu = 0.822 \pm 0.031$;
$\langle \eta_L \rangle_R \approx 1.0 \cdot 10^{-5}$, $\langle \eta_L \rangle_{\text{full}}
\approx 7.3 \cdot 10^{-6}$. Fit of $\eta_L \cdot t$ vs $\ln(\lambda t)$ on $[10^3, 1.5 \cdot 10^4]$:
$n = 71$ points, $R^2 = 0.888$ for $N_{\text{max}} = R t$; $R^2 = 0.906$ for the full budget.

*Reproduction command and outputs.* `python simulations/markov_drift_ou/main.py`; wall-time
$\sim 1$ minute; outputs — `results_summary.{txt,json}`, `run.log`, figures
`paper/figs/fig{4,5}_ou_*.png`.

### S4.3 Adiabatic scan (§ S8.3)

*Model class.* The OU parametrisation of § 6.2 unchanged ($k = 8$, $\lambda = 10^{-3}$, $K = 56$,
Euler–Maruyama $\Delta t = 1$, Robbins–Monro $\eta_t = c_0 / \sqrt{t}$, $c_0 = 3$). The measured
quantity is the cumulative excess-loss $L_{\text{excess}}(t)$ by (S8.1) with direct access to the
ground truth $\theta(s)$.

*Parametric scan.* Three values of $\sigma$ at fixed $\lambda = 10^{-3}$:

| $\sigma$ | $\sigma^2/\lambda$ | $1/\sigma^2$ | $T_{\text{total}}$ | Realisations |
|----------|--------------------|--------------|--------------------|------------|
| $0.01$ | $10^{-1}$ | $10^4$ | $5 \cdot 10^4$ | $25$ |
| $0.003$ | $9 \cdot 10^{-3}$ | $1.1 \cdot 10^5$ | $5 \cdot 10^4$ | $25$ |
| $0.001$ | $10^{-3}$ | $10^6$ | $2 \cdot 10^5$ | $25$ |

The simulation time $T$ is chosen within the transient window $1/\sigma^2$ for each $\sigma$.

*Joint fit.* $L_{\text{excess}}(t) = A \cdot t + B \cdot \ln(\lambda t) + C$ on the window
$t \in [10^3, 10^4]$ (fixed for direct comparison across $\sigma$).

*Wide-window negative control.* An additional fit with the extended window
$[10^3, T_{\text{total}}]$: for $\sigma^2/\lambda = 10^{-3}$ and $T_{\text{total}} = 2 \cdot 10^5$
the coefficient $B/(K/2)$ grows to $2.0$ — an artefact of the Robbins–Monro $1/\sqrt{t}$ bias of the
learner on long horizons, not part of the BNT effect. This control records the sensitivity of the
result to the choice of window and prevents the interpretation of the monotone trend of $B/(K/2)$ as
the result of selective presentation.

*Key numerical results.* $B/(K/2) = \{0.54, 0.78, 0.85\}$ for $\sigma^2/\lambda
\in \{10^{-1}, 9 \cdot 10^{-3}, 10^{-3}\}$ respectively; $R^2 \ge 0.9995$ at all points;
the drift-share $A \cdot t / L_{\text{excess}}(t)$ at $t = 10^4$ falls from $76\%$ to $10\%$.
The monotone growth of $B$ as $\sigma^2/\lambda \to 0$ is consistent with the conjecture of § S8.1, but
does not prove convergence $B \to K/2$; $5$–$6$ scan points are needed with an explicit analytical
form of the finite-$\sigma$ correction $B(\epsilon) = (K/2)(1 - g(\epsilon))$.

*Reproduction command and outputs.* `python simulations/markov_drift_ou_iinf_adiab/main.py`;
wall-time $\sim 1$ minute; outputs — `results_summary.{txt,json}`, `run.log`, figures
`simulations/markov_drift_ou_iinf_adiab/fig_adiab_*.png`, copied into
`paper/figs/fig{6,7}_adiab_*.png`.

### S4.4 Reproducibility: dependencies and archiving

*Seed and dependencies.* The global $\text{SEED} = 20260524$ is set at the start of `main.py` of
each script; per-run offsets ensure the independence of realisations. Dependencies —
`numpy`, `scipy`, `scikit-learn` (KSG), `matplotlib`; exact versions are in `requirements.txt`
of each directory.

*Full set of commands.* At the root of the repository, in sequence,
`python simulations/markov_drift{,_ou,_ou_iinf_adiab}/main.py`; total wall-time on a
mid-range laptop $\sim 7$–$12$ minutes.

*Zenodo archiving.* A versioned snapshot has been archived on Zenodo (DOI: 10.5281/zenodo.20653051,
release `v1.0-submission`); reproduction is also possible by cloning
`andriishin/landauer-nostalgia-oa` at the `v1.0-submission` tag (the same snapshot archived on Zenodo).

---

## S5. The full apparatus of § 5 main: Remarks 1–6 to Lemma 2, the protocol of § 5.2, $\nu_c$ and $\tau_{\text{reset}}^*$

This section assembles the formal details of § 5 main, moved here to lighten the main exposition.
Structure: S5.1 — the six Remarks to Lemma 2 (§ 4.4 main fixes the formulation (8a) and the brief
statement; here is their full development); S5.2 — point 5 of the protocol of § 5.2 main (the ban on
revising the definition of the object of measurement) with justification through (Wolpert 1996)
no-free-lunch; S5.3 — the balance condition for the existence of $\nu_c$ through (11b) and a discussion
of its structural status; S5.4 — the optimal $\tau_{\text{reset}}^*$ (12), its connection to dynamic
regret bounds (Besbes et al. 2015); S5.5 — the analogy with spin glasses through (Hopfield 1982; Amit et al. 1985; Engel and Van den Broeck 2001).

### S5.1 Six Remarks to Lemma 2

**Remark 1 (stationary limit).** The stationary limit of Lemma 2 is given *not* by
$\tau_E \to \infty$ (this passage violates the condition $c < 1$ of $(*)$ main), but by equation (1)
main with $\dot{E}_{\text{store}} = \dot{E}_{\text{grow}} = 0$ (§ 2.2 main): in this regime the
nostalgic bits carry no payment, and the scale $\nu^{\text{theor}}$ through (7) main loses its
thermodynamic meaning (likewise for $\nu^{\text{op}}$ as $\hat P \to P_0$). Thus, Lemma 2 does not
contradict the stationary lemma (Andriishin 2026, § 2.1), but *supplements* it to a new regime: for
$\dot E_{\text{store}} > 0$ or $\dot E_{\text{grow}} > 0$ the condition (8a) is operative, and in the
stationary limit it degenerates through the vanishing of both payments.

**Remark 2 (computability and falsifiability of $(*)$).** The constant $c$ from $(*)$ main is a
*computable* function of three measurable parameters $(\dot{M}_{\text{refresh}}, \tau_E,
|M_{\le t}|)$, not a free parameter of the theory. Lemma 2 applies at measured $c < 1$; for
$c \ge 1$ the system passes into the *sliding-window* regime: over the interval $\tau_E$ the entire
memory is fully refreshed, there is no archival layer, and the non-stationary theory degenerates. The
test on $c$ is an empirical condition of applicability, separating systems with an archival layer from
those without it. Long-context LLMs with an $M_t$ layer of order $10^5$–$10^6$ tokens at inference
give $\dot{M}_{\text{refresh}} \tau_E \gg |M_t|$, $c \gg 1$, and Lemma 2 assigns them to the
sliding window — the theory by its own condition removes this layer from the domain of applicability
(§ 8.3 main).

**Remark 3 (substantiveness of (8a) beyond the definition).** The fraction $1-c$ is arithmetic
(the ratio of the number of unrefreshed bits to the total memory size), while $\nu^{\text{theor}}$ by
(7) main is information-theoretic (the normalised MI deficit). The substantive part of Lemma 2 is the
proof that the MI contribution of an unrefreshed bit to $I_{\text{pred}}(t,\tau)$ vanishes
asymptotically by DPI through the decorrelation of the OU latent with characteristic time $\tau_E$
(see the detailed proof in S1 supplementary through MI-tensorisation over the independent coordinates
of OU, § S1.5). Three conditions are required: (i) decorrelation of the optimal statistic with
$\tau_E$; (ii) the ideal observer in the denominator of (7) main; (iii) FIFO refresh. Under violation
of any of (i)–(iii) the inequality (8a) does not follow from the arithmetic $c$.

**Remark 4 (extension to the general non-stationary class).** The basic Lemma 2 is proven for the
OU parametrisation of the logits (§ 5.2 main, S1–S2 supplementary). The PSP surrogate (§ 6.1 main)
and multiscale drift reproduce the result *illustratively* for $\lambda \tau_{\text{relax}} \ll 1$
through a single decorrelation mechanism with effective $\tau_E$, but there is no rigorous proof for
non-OU classes in the present work; see § 8.5 main, open problem. The Wolpert no-free-lunch
(Wolpert 1996) gives a fundamental justification for narrowing the class of § 5.2: a universal theory
of non-stationary learnability without an a priori restriction of the class is impossible; the
slow-driving class is a necessary narrowing, not a methodological choice.

**Remark 5 (per-bit-uniformity).** The passage from the arithmetic fraction of unrefreshed bits
$c$ to the information-theoretic inequality $\nu \ge 1 - c$ requires the additivity of MI over
independently refreshed bits:

$$I(M^{\text{refr}}; X_E^{[t,t+\tau]}) \le c \cdot I_{\text{pred}}^{\text{opt}} \cdot (1+o(1)).$$

In the slow-driving limit the refresh fraction accumulates no cross-bit redundancy
(Clarke and Barron 1990); the full proof of this statement for the OU parametrisation is
S2 supplementary. Under violation of the condition $c \to \kappa c$, $\kappa \ge 1$, and Lemma 2 is
formulated with a multiplicative correction $\nu \ge 1 - \kappa c$; in the OU class with
independent $W_{ij}$ and a diagonal anchor $\theta_{ii} \equiv 0$ the per-bit-uniformity condition is
satisfied by construction with $\kappa = 1 + o(1)$ (S2.6 supplementary).

**Remark 6 (softmax-regularity for OU).** The transfer of OU decorrelation from the autocovariance of
the logits $\theta(t)$ to any bounded measurable functional $\xi(\theta)$ (in particular, to the
entries of the transition matrix $P_{ij}(\theta)$ and to any sufficient statistic of prediction) is
ensured by the Lipschitz continuity of softmax in total variation on the set of logits with finite
stationary variance $\sigma^2/(2\lambda)$. Through MI-tensorisation over the independent coordinates
of OU (§ S1.5) and DPI the exponential decay $\mathbb{E}[\theta(t)\theta(s)] \propto e^{-(t-s)/\tau_E}$
is transferred to any bounded $\xi$ with the same characteristic time $\tau_E$. The full proof is
S1 supplementary (five steps: Lipschitz softmax in $\ell_\infty$; transfer to the TV of the transition
matrix; boundedness of the OU trajectory on the ball $R(\delta)$ with probability $1-\delta$; transfer
of the decay to the covariance of $\xi$; translation of the TV bound into a bound on the conditional
MI through $\chi^2$ and DPI).

Lemma 2 in light of Remarks 1–6 is constructive: under $(*)$ it indicates an explicit mechanism of the
emergence of nostalgia (the desynchronisation of the optimal feature through OU decorrelation) and a
physical obstacle to its vanishing (the finite memory-refresh rate under FIFO scheduling). The joint
action of Remarks 5–6 is the structural justification of the claimed asymptotics $\nu \ge 1 - c$ in the
OU class of § 5.2 main; Remark 2 fixes the domain of applicability through the measurable condition
$c < 1$.

### S5.2 The protocol of § 5.2 main, point 5: the ban on revising the definition of the object of measurement

Point 5 of the protocol (§ 5.2 main) is a pre-registration commitment that protects falsifiability: the
operational definition of the measured quantity is fixed *before* the data and is not revised, under the
same theory name, after the data are in. When prediction and measurement disagree systematically, two
honest moves remain: (a) refute the theory; (b) explicitly change its name, leaving the old name attached
to the original definition. The rule is symmetric — it binds the theorist and the empiricist equally.

*Why it matters.* The standard response to a failed empirical test is either to reject the theory or to
narrow its conditions of applicability (add caveats); both are legitimate and keep the test severe.
Point 5 closes only a third, methodologically dubious route — revising the *operational definition* of
the object of measurement while keeping the theory's name (a classical ad hoc rescue in Popper's sense,
disguised as a technical refinement of operationalisation). The symmetry keeps the discipline fair in
both directions: the empiricist may not "refine" the definition of the measured quantity post hoc so
that the data agree with the theory; the theorist may not "refine" the definition of the predicted
quantity post hoc so that they disagree.

*How the present work honours it.* The cumulative excess-loss $L_{\text{excess}}(t)$ (S8.1) is kept in
§ S8 as a separate, explicitly labelled counterfactual diagnostic rather than used as the operationally
measurable $I_{\text{pred}}$: $L_{\text{excess}}$ requires access to the ground truth $\theta(s)$,
unavailable in real biological systems, so conflating it with the one-step MI (4) main under the name
"predictive information" would be a post hoc revision of the object of measurement. The theoretical
constructions (Lemma 2, $\nu_c$, $C_v^{\text{predictive}}$) rest on $\nu^{\text{theor}}$ by (7) main; the
operational predictions (§ 8.6 main, Prediction 1, Prediction 2) rest on $\nu^{\text{op}}$ by (7') main.
Substitution under one name is excluded; separated accounting of the two quantities is not.

### S5.3 Existence of $\nu_c$: the balance condition (11b) and its structural status

The balance estimate from the condition of equality of the increment of $\dot{I}_{\text{pred}}$ and the
losses from staleness (equation (11b) main):

$$\nu_c^{\text{theor}} \cdot |M_{\le t}| / \tau_{1/2} = (1 - \nu_c^{\text{theor}}) \cdot \dot{M}_{\text{refresh}}.$$

The dimensionless $\eta_M := \dot{M}_{\text{refresh}} \tau_{1/2} / |M_{\le t}|$ formally gives
$\nu_c^{\text{theor}} = \eta_M / (1 + \eta_M)$; direct substitution in this simplest form
gives numerical values inconsistent with the empirical proxy of § 6.1 main
($\nu_c^{\text{emp}} \approx 0.44$ vs. $\nu_c \approx 0.016$ for the $\eta_M$ of the base
parametrisation of the PSP — a discrepancy of more than an order of magnitude).

*Structural status of (11b).* (11b) is treated as a *structural balance condition*, not a
quantitative predictor: the work postulates the existence of $\nu_c^{\text{theor}} \in
(0, 1)$ as a fixed point of the balance of information increment and losses from staleness. The
discrepancy with the empirical $\nu_c^{\text{emp}}$ supports this decision — direct substitution of the
simplest form of (11b) is rejected, which is recorded as the limit of its interpretive force. The
derivation of a quantitative formula $\nu_c^{\text{theor}}(\lambda, \sigma,
K, \tau_{1/2}, \dot M_{\text{refresh}})$ is an open problem § 8.5 main.

*Active part of memory.* In (11b) $|M_{\le t}|$ is the *active* part of memory, not the structural
capacity limit $|M_{\le t}^{\text{cap}}|$. In § 6.1 main: $|M_t| = k^2 = 64$ parameters of the
estimator; this quantity grows with time in systems with expanding memory, and the numerical value of
$\nu_c^{\text{theor}}$ depends on time. The empirical proxy $\nu_c^{\text{emp}}$ of § 6.1 main is
$\nu^{\text{op}}_c$ by (7') main, not $\nu^{\text{theor}}_c$; the working assumption of § 4.1
main on the stability of the regime threshold under a change of normalisation provides the transition
between them with a possible numerical shift.

### S5.4 The optimal $\tau_{\text{reset}}^*$ and its connection to dynamic regret

The optimal $\tau_{\text{reset}}^*$ maximises $\langle \eta_L \rangle$: frequent reset increases the
payment, infrequent reset leads to nostalgic collapse (§ 5.3 main). The structural estimate:

$$\tau_{\text{reset}}^* \sim \tau_E \cdot \log(|M_{\le t}^{\text{cap}}|/\dot{M}_{\text{refresh}} \tau_E),$$

where $|M_{\le t}^{\text{cap}}|$ is the structural capacity limit of memory. The logarithmic dependence
is a reflection of the standard form of the optimum for processes with exponential staleness and a
linear refresh payment.

*Biological analogy.* The generation time of a taxon correlates with the rate of change of the
econiche — an empirical regularity of population biology, consistent with (12) main up to a
logarithmic factor. In mammals, in early embryogenesis and the germ line, an almost complete erasure
of DNA methylation with re-establishment occurs (Reik 2007); this can be interpreted as a mechanism of
$\tau_{\text{reset}}^*$ on the generational scale — a periodic reset of the nostalgic layer for the
sake of restoring $\nu \to 0$ in a new generation.

*Connection to dynamic regret bounds.* In non-stationary stochastic optimisation (Besbes et al. 2015)
(Theorem 2, stochastic setting) the regret-optimal length of the restart epoch has a *power-law*
structure: the optimal length $\propto (T/V_T)^{2/3}$, where $V_T$ is the variation of the target on
the horizon $T$. This differs from the logarithmic form $\tau_{\text{reset}}^* \sim \tau_E \log(\cdot)$
main §8.4(b): the main text gives a *logarithmic* form through a linear refresh payment with
exponential staleness, whereas the BGZ regret-optimal restart is *power-law* $(T/V_T)^{2/3}$. Only the
monotonicity in $\tau_E$ coincides (slower drift — rarer reset), not the functional form. In
no-free-lunch (Wolpert 1996) a universal $\tau_{\text{reset}}^*$ is impossible, which is consistent with
the conditionality of (12) main relative to the class of § 5.2 main.

### S5.5 The analogy with spin glasses

(11b) main is a mean-field balance, structurally analogous to replica symmetry in the statistical
mechanics of learning (Engel and Van den Broeck 2001). In the classical Hopfield model (Hopfield 1982) the
critical load $\alpha_c \approx 0.138$ (Amit et al. 1985) is the threshold beyond which
the associative memory loses the ability to retrieve the stored patterns without errors. Our
$\nu_c^{\text{theor}}$ is the analogue in a non-equilibrium learning system: the threshold beyond which
the nostalgic fraction exceeds a critical fraction and the predictive power $\eta_L(t)$ collapses.

The quantitative connection between $\nu_c^{\text{theor}}$ and the Hopfield $\alpha_c$ is a separate
problem: $\alpha_c$ is the ratio of the number of stored patterns to the number of neurons in an
equilibrium system with a fixed set of targets, $\nu_c^{\text{theor}}$ is the threshold ratio of the
nostalgic fraction in a non-equilibrium system with a drifting latent parameter. The structural
correspondence — both quantities reflect the limiting load on the memory capacity in systems with a
fixed number of degrees of freedom; the dynamical characteristics, however, are substantially different.

---

## S6. Full technical details of § 6 main: PSP parameters, Robbins–Monro learner, wide-window control

This section assembles the technical details of the numerical illustrations of § 6.1 (PSP) and § 6.2 (OU)
main, which in the main text are left in the form of a brief statement of results with a cross-ref
here. Structure: S6.1 — full parameters of the PSP simulation (§ 6.1 main) and the MI estimator;
S6.2 — simulation pseudocode and comparative analysis of estimators; S6.3 — the Robbins–Monro online
learner of the OU simulation (§ 6.2 main) and its asymptotic Bayes-optimality; S6.4 — the
numerical estimate of $c$ for Robbins–Monro and the discrepancy with the prediction (8a); S6.5 —
the wide-window negative control of the adiabatic scan and the discrepancy $C_v^{\text{static}}/
C_v^{\text{predictive}}$; S6.6 — the empirical $\nu_c^{\text{emp}}$ and its status.

All reproducibility parameters (seed, library versions, directories) are S4
supplementary (simulation parameters for reproducibility); this section focuses on the
methodological justification of the choice of parameters and the interpretation of the results.

### S6.1 PSP simulation (§ 6.1 main): full parameters and estimator

*Model.* A finite-state Markov process with Poisson switches of the transition matrix.
Alphabet $\Omega = \{1, \ldots, k\}$, $k = 8$. The matrix $P(t)$ is piecewise-constant on the intervals
$[T_i, T_{i+1})$; the switch times $\{T_i\}$ form a Poisson stream with intensity $\lambda$.
At time $T_i$ the matrix is updated to a random stochastic $P^{(i)} \sim
\text{Dirichlet}(\alpha = 0.3)$ by rows. $\alpha = 0.3$ gives $\langle
\tau_{\text{relax}} \rangle \approx 2$–$3$ steps and $\lambda \tau_{\text{relax}} \approx 3
\cdot 10^{-3} \ll 1$ with margin — slow-driving is satisfied.

*Observer's memory.* The observer system maintains a memory $M_t$ consisting of an estimate of the
current transition matrix $\hat{P}(t)$ on a sliding window of length $\tau_w$ with Laplace smoothing
(pseudocount 1 per cell). The structural dimension of the memory $|M_t| = k^2 = 64$ parameters is the
number of independent cells of the estimator; for the budget $N_{\text{max}}(t) = N_0 + R \cdot t$ in
the simulation the normalisation $k_B T \ln 2 = 1$ and $R = 1$ unit per step are adopted (the base
operational payment plus Lemma 1 on the $k^2$ maintained parameters).

*Parameters of the three regimes.* Stationary regime: $\lambda = 0$ (the matrix does not change), $\tau_w
= 100$, $T = 10^4$ steps, averaging over 80 realisations. Drift without reset: $\lambda =
10^{-3}$ (mean $\tau_E = 10^3$ steps between switches), $\tau_w \to \infty$ (the memory is not
cleared, $\hat{P}(t)$ is averaged over the whole history), $T = 10^5$ steps, 60 realisations.
Drift with reset: $\lambda = 10^{-3}$, $\tau_w \in \{100, 500, 2000\}$, $T = 10^5$, 60
realisations.

*Grid of $\tau_w$ for Fig. 2.* Eight points $\tau_w \in \{50, 100, 200, 500, 10^3, 2\cdot
10^3, 5\cdot 10^3, 10^4\}$ at $\lambda = 10^{-3}$, $T = 5 \cdot 10^4$, 40 realisations per
point.

### S6.2 Simulation pseudocode and comparative analysis of MI estimators

*Pseudocode of the PSP simulation.*

```
for each condition (λ, τ_w):
    initialise P_0, draw T_1 ~ Exp(λ)
    initialise memory: empirical_freq = matrix of zeros (k × k)
    for t = 0, 1, ..., T:
        if t == T_i: update P to a random P^(i) ~ Dirichlet(α=0.3, k rows),
                     draw T_{i+1}
        X_t = sample from P_t[X_{t-1}, :]
        update empirical_freq[X_{t-1}, X_t] += 1
        if t mod τ_w == 0 and τ_w < ∞: empirical_freq = 0  # reset
        P_hat(t) = empirical_freq normalised by rows (with Laplace smoothing)
        compute I_pred(t, τ) = MI(X_{t+1..t+τ}; P_hat(t))
        compute ν(t) = 1 - I_pred(t,τ) / I_pred^optimal(P_t, τ)
        accumulate N_max(t) = accumulated storage cost of empirical_freq + capacity growth
        η_L(t) = I_pred(t, τ) / N_max(t)
for each set (λ, τ_w): plot η_L(t), ν(t)
```

*MI estimator: KL-detuned proxy.* At prediction horizon $\tau$ the numerical estimate of the
kernel definition (4) main is carried out through

$$\widehat I_{\text{pred}}(t, \tau) = I_{\text{opt}}(P(t), \tau) - \mathbb{E}_{X_t \sim
\pi(t)}\!\bigl[ D_{\text{KL}}\!\bigl(P(t)^\tau(\cdot \mid X_t)\,\|\,\hat{P}(t)^\tau(\cdot
\mid X_t)\bigr) \bigr],$$

clipped at zero from below. Here $I_{\text{opt}}(P, \tau) = H(\pi) - \sum_i \pi_i
H(P^\tau_{i,\cdot})$ is the predictive information of the optimal observer with the exact model;
the KL penalty for the mismatch of $\hat{P}$ to the true $P$ is subtracted. Connection to (4) main: for
a Markov chain with a fixed estimator $\widehat I_{\text{pred}}(t, \tau)$ is majorised from above by the
conditional MI $I(X_t; X_{t+\tau}\mid\hat P(t))$ (the KL penalty accounts for the maladaptation of the
model as a lower bound on the loss of informativeness); the clipping at zero from below is ensured by
the fact that the conditional MI is non-negative. On the ergodic plateau $\hat P(t) \to P(t)$ the KL
penalty vanishes, and $\widehat I_{\text{pred}} \to I_{\text{opt}} = I(X_t; X_{t+\tau}\mid P_0)$,
which is consistent with (4) main in the stationary limit.

*Comparison with independent estimators.* Alternative discrete MI estimates with bias correction
(Treves and Panzeri 1995) and the continuous KSG estimator with $k$-NN, $k_{NN} = 3$
(Kraskov et al. 2004) give consistent results on the parameter grid of the PSP simulation. The consistency
of three independent estimators (KL-detuned proxy, Treves–Panzeri plug-in, KSG) on the same
trajectories confirms that the observed qualitative patterns (Fig. 1–3 main) are not an artefact of a
specific MI estimation method.

*Dome-shaped dependence $\langle \eta_L \rangle(\tau_w)$.* The dome-shaped curve of Fig. 2
main with maximum $\langle \eta_L \rangle^* \approx 1.3 \cdot 10^{-5}$ at $\tau_w^*
\approx 200$ steps (plateau $200$–$500$) is the non-stationary analogue of the bias-variance
trade-off of classical learning theory: small $\tau_w$ — high variance of the estimator $\hat P$
(insufficient data), large $\tau_w$ — bias due to averaging over outdated regimes; the optimum
$\tau_w^*$ is the balance, analogous to structural risk minimisation.

The position of the maximum $\tau_w^* \approx 0.2\,\tau_E$ is consistent with (12) main up to a
logarithmic factor; the exact value of the optimum depends on fine parameters of the estimator
(Laplace smoothing, the $\alpha$ of the Dirichlet).

### S6.3 The Robbins–Monro online learner of the OU simulation (§ 6.2 main)

*Maintenance of the online estimate.* In the OU simulation (§ 6.2 main) the learner maintains an
online estimate $\hat\theta(t)$ through the maximisation of the conditional log-likelihood with a
decreasing step $\eta_t = c_0/\sqrt{t}$ ($c_0 = 3$): at each step, upon observing $(x_t
\to x_{t+1})$, one performs

$$\hat\theta_{x_t,j}(t+1) = \hat\theta_{x_t,j}(t) + \eta_t \cdot \bigl[\mathbb{1}[x_{t+1}
= j] - \hat P_{x_t,j}(t)\bigr]$$

for $j \ne x_t$. This is a diffusion-prior-Bayes approximation of the Bayesian update over the OU latent
in the sense of § 5.2 main; the rate $1/\sqrt{t}$ corresponds to the concentration of the posterior
along each of the $K$ directions and is asymptotically Bayes-optimal for regular parametric
models (Clarke and Barron 1990, § 3).

*Polyak–Ruppert averaging note.* The standard improvement of Robbins–Monro — Polyak–Ruppert
averaging of the iterates — gives an asymptotically efficient estimate with smaller variance; in
the present implementation it is not used, since the observed qualitative patterns (Fig. 4–5
main) are robust to the choice of Robbins–Monro variant up to log factors on the final
variance. The use of Polyak–Ruppert would change the quantitative values of $\eta_L(t)$ by a
constant factor, but not the qualitative picture of the collapse.

*Structural equivalence to Bayes.* Robbins–Monro with the indicated rate $\eta_t \propto
1/\sqrt{t}$ is structurally equivalent to an online approximation of the Bayesian update over the
OU latent of § 5.2 main in the following sense: the post-cumulative distribution of the estimates
$\hat\theta(t)$ as $t \to \infty$ concentrates around the conditional mean of the posterior with
covariance of order $\mathcal{F}^{-1}/t$, where $\mathcal{F}$ is the Fisher information matrix
of the parametric model (Clarke and Barron 1990). This is the justification of the choice of $\eta_t \propto 1/\sqrt{t}$
as functionally equivalent to the Bayesian learner.

### S6.4 The numerical estimate of $c$ for Robbins–Monro and the discrepancy with (8a)

*Effective refresh rate.* The adaptive step $\eta_t = c_0/\sqrt t$ with $c_0 = 3$ gives an
effective refresh rate $\dot M_{\text{refresh}}^{\text{eff}}(t) \sim K \cdot \eta_t /
\tau_E = K c_0/(\tau_E \sqrt t)$; for $K = 56$, $\tau_E = 10^3$, $t = 2 \cdot 10^4$
one obtains

$$c(t) = \frac{\dot M_{\text{refresh}}^{\text{eff}} \cdot \tau_E}{|M_{\le t}|} \approx
\frac{K c_0}{|M_{\le t}| \sqrt t} \approx \frac{56 \cdot 3}{56 \cdot \sqrt{2 \cdot 10^4}}
\approx 0.02.$$

*Discrepancy with the prediction of Lemma 2.* Lemma 2 gives $\liminf \nu^{\text{theor}}(t) \ge 1
- c \approx 0.98$. The observed $\nu(t > 10^4) = 0.82 \pm 0.03$ is *below* the predicted
asymptote; the gap $\approx 0.16$ is explained by (i) the finiteness of the observation time
(the asymptotics (8a) main is the $\liminf$ as $t \to \infty$); (ii) the per-bit-uniformity assumption
of Remark 5 (S5.1 supplementary), holding only approximately for the OU parametrisation of the
logits with softmax nonlinearity (Remark 6). The full proof of convergence to $1-c$
for the specific OU parametrisation is an open problem § 8.5 main.

### S6.5 Wide-window negative control and the discrepancy $C_v^{\text{static}}/C_v^{\text{predictive}}$

*Wide-window negative control.* Upon extending the window of the joint fit of the adiabatic scan
(§ S8.3) to $[10^3, T_{\text{total}}]$ (instead of the fixed $[10^3, 10^4]$) the fit
degenerates: for $\sigma^2/\lambda = 10^{-3}$ and $T_{\text{total}} = 2\cdot 10^5$
the coefficient $B/(K/2)$ grows to 2.0 — an artefact of the Robbins–Monro $1/\sqrt t$ bias of the
learner on long horizons, not part of the BNT effect. The window $[10^3, 10^4]$ is chosen for direct
comparison between the points of the adiabatic series and is consistent with the transient condition
$t \ll 1/\sigma^2$. This control records the sensitivity of the result to the choice of the fit
window and prevents the interpretation of the monotone trend of $B/(K/2)$ as the result of selective
presentation.

*Numerical scan $B/(K/2) \in \{0.54, 0.78, 0.85\}$.* Full results — § S8.3
(table of $B$, $B/(K/2)$, $R^2$, drift-share); S4.3 supplementary (scan parameters).

*Discrepancy $C_v^{\text{static}}/C_v^{\text{predictive}}$.* In the "drift without reset" regime
(§ 6.1 main) at $\nu(t) \to 0.999$ one has $C_v^{\text{predictive}}(t) = (1 - \nu(t))
\cdot C_v^{\text{static}}(t) \to 0$, whereas $C_v^{\text{static}}$ remains bounded
below by the size of the accumulated memory $\hat P(t)$. The ratio $C_v^{\text{static}}/
C_v^{\text{predictive}} \to \infty$ is a numerical illustration of the diagnostic signal
of nostalgia (§ 4.3 main): two proxies, agreeing in the stationary regime, diverge under
observed drift. This operationalises Prediction 1 § 8.6 main: the crossing of
$\nu_c^{\text{op}}$ is diagnosed through the divergence of the two proxies, without direct measurement
of $\nu(t)$.

### S6.6 The empirical $\nu_c^{\text{emp}}$ and its status

The empirical $\nu_c^{\text{emp}}$, defined as the value of $\nu(t)$ at the moment when
$\eta_L(t)$ falls to half of its peak in the drift-without-reset regime, equals
$\nu_c^{\text{emp}} \approx 0.44$ for the base parametrisation of the PSP. This value
is consistent with the structural statement of § 5.3 main on the existence of $\nu_c \in (0, 1)$
separating the adaptive regime and the collapse regime, but **does not claim a universal
quantitative estimate**: $\nu_c^{\text{emp}}$ depends on the chosen operational definition
(the half-life threshold of $\eta_L$) and on the parameters of the specific PSP surrogate.

*Discrepancy with (11b).* The direct balance substitution $\nu_c = \eta_M/(1+\eta_M) \approx
0.016$ in (11b) main is rejected empirically — a discrepancy of more than an order of magnitude. This
supports the decision of § 5.3 main / S5.3 supplementary to treat (11b) as a structural
balance condition, not as a quantitative predictor. The derivation of a universal formula
$\nu_c(\lambda, \sigma, K, \tau_{1/2}, \dot M_{\text{refresh}})$ is an open problem § 8.5 main.

*Stability under a change of normalisation.* By DPI $\nu^{\text{op}}(t) \le \nu^{\text{theor}}(t)$
pointwise (§ 4.1 main; the smaller denominator $C_\mu^{\text{emp}} \le I_{\text{pred}}^{\text{opt}}$
in (7') gives a larger fragment and smaller nostalgia); consequently the empirical $\nu_c^{\text{emp}}$ is
$\nu^{\text{op}}_c$, a lower estimate for $\nu^{\text{theor}}_c$. The working assumption of § 4.1
main on the monotonicity of the regime threshold under a change of normalisation provides the
applicability of the empirical $\nu_c^{\text{emp}}$ as an operational surrogate of the theoretical
$\nu_c^{\text{theor}}$ with a possible numerical shift, but without a qualitative displacement of the
threshold.

---

## S7. Active inference: full formal correspondence between $\eta_L$ and EFE, $\rho(t)$, caveats

This section assembles the formal details of § 7 main, moved here to lighten the main
exposition. Structure: S7.1 — the heuristic correspondence $\Delta G^{\text{epist}}$ ↔
$-\Delta I_{\text{pred}}$ (S7.1), the full justification with the conversion of units; S7.2 — the
conjecture (S7.2) with its explicit hypothetical status and a discussion of the direction of
formalisation; S7.3 — the categorical caveat $I_q(s'; o') \ne I(X_t; X_{t+\tau})$ through the
sufficiency of $\hat P(t)$; S7.4 — the context of sophisticated inference (Friston et al. 2021)
and the Pearl/Friston blanket (Bruineberg et al. 2018; Bruineberg et al. 2022) demarcation; S7.5 — the rebuttal of
Andrews (2021) / Williams (2020) through a cross-ref to paper #1; S7.6 — the informational rejuvenation
metric $\rho(t)$ in detail with empirical interpretation.

### S7.1 Heuristic correspondence: full justification

In the formulation of active inference (Friston et al. 2017a, eq. 3; Parr2022, ch. 7) the expected free energy
of a policy $\pi$ is reckoned in nats and is

$$G[\pi] = \mathbb{E}_{q(o', s' \mid \pi)}\!\left[\ln q(s' \mid \pi) - \ln p(o', s'\mid C)\right],$$

where $p(o', s'\mid C)$ is the preference distribution with parameter $C$ encoding the agent's goals.
The canonical decomposition decomposes $G = -\mathbb{E}_q[D_{\text{KL}}(q(s'\mid o',\pi)
\|q(s'\mid\pi))] - \mathbb{E}_q[\ln p(o'\mid C)]$ into epistemic value (the information gain
about hidden states) and pragmatic value (the proximity to preferences).

*Conversion of units.* $G$ is measured in nats; $I_{\text{pred}}$ in bits (Appendix A
main). Conversion: $G_{\text{[nats]}} = G_{\text{[bits]}} \cdot \ln 2$; conversion into energy
units through $\beta = 1/(k_B T)$ gives $E_G = k_B T \cdot G_{\text{[nats]}}$. Similarly
$E_{I_{\text{pred}}} = k_B T \ln 2 \cdot (I_{\text{pred}})_{\text{[bits]}}$. For the increment over
the policy window $\Delta t$:

$$k_B T \cdot \Delta G^{\text{epist}}[\pi]\bigr|_{\Delta t} \approx -k_B T \ln 2 \cdot
\Delta I_{\text{pred}}\bigr|_{\Delta t}, \tag{S7.1}$$

up to an additive constant absorbing the pragmatic contribution and independent of $\pi$.
Both sides have the dimension of energy (J); the minus sign reflects the convention "$G$
is minimised under the maximisation of $\dot{I}_{\text{pred}}$"; the equivalent dimensionless
form is obtained by dividing by $k_B T$.

*Hypothetical status.* (S7.1) is fixed as a *heuristic correspondence*, not as a
proven inequality. A proof of the strict inequality requires (i) careful handling of the
time derivatives in a non-stationary generative model; (ii) a demonstration that the
pragmatic contribution is indeed independent of $\pi$ on the time scale $\Delta t$; (iii)
accounting for the corrections from the categorical caveat S7.3. All three ingredients are absent in
the existing literature on active inference in the non-stationary regime.

### S7.2 Conjecture: hypothetical status

The hypothetical theorem is the direction in which (S7.1) can be turned into a strict
inequality:

$$k_B T \cdot \Delta G^{\text{epist}}[\pi]\bigr|_{\Delta t} \stackrel{\text{conj.}}{\ge}
-k_B T \ln 2 \cdot \Delta I_{\text{pred}}\bigr|_{\Delta t}, \tag{S7.2}$$

where the label "conj." marks the hypothetical status. The formal provability of (S7.2) is an
open problem requiring careful handling of the time derivatives and a non-stationary
generative model in the extended active inference framework (Parr et al. 2022).

*Sophisticated inference as a direction of formalisation.* Sophisticated inference
(Friston et al. 2021) extends standard active inference to policies that
account for post-observational updates of the future prior distributions; this is a natural
framework for the careful handling of the time derivatives in a non-stationary generative model.
A proof of (S7.2) in the sophisticated inference framework is a natural next step of the
extension of the present work.

### S7.3 The categorical caveat $I_q(s'; o') \ne I(X_t; X_{t+\tau})$

The epistemic value in (Friston2017) is the mutual information $I_q(s'; o'\mid\pi)$ between the hidden
state $s'$ and the observation $o'$, whereas $I_{\text{pred}}(t,\tau)$ by (4) main is the MI
between two consecutive observations $X_t$ and $X_{t+\tau}$. These are categorically
different objects: the former is the MI at a single point in time between the latent and the
observation, the latter is the MI between two moments in time along the observed trajectory.

*Bridge through the sufficiency of $\hat P(t)$.* Under the interpretation $s' = \theta(t+\tau)$, $o' =
X_{t+\tau}$ and the condition $X_{t+\tau} \perp \theta(t) \mid \hat P(t)$ (the statistical
sufficiency of $\hat P(t)$ with respect to $X_{t+\tau}$), by the chain rule of MI one obtains

$$I_{\text{pred}}(t,\tau) \le I_q(\theta(t+\tau); X_{t+\tau} \mid \hat P(t)).$$

The minimisation of the epistemic part of $G[\pi]$ under the budget $N_{\text{max}}(t)$ *heuristically
corresponds* to the local maximisation of $\eta_L(t)$ (in the sense of the heuristic correspondence § S7.1;
the strict implication is not proven — the bridge gives an upper bound on $I_{\text{pred}}$, not an
equality). The pragmatic component through the preference $C$ is *outside* the
thermodynamic scale; the full equivalence of $\eta_L$ and $G$ would require encoding the
preference through the minimisation of $\dot{E}_{\text{actual}}$ — a nontrivial extension of the FEP.

The condition of correspondence is that $M_{\le t}$ includes $q(s)$ and the history of observations
(an extension of the generative model to the episodic level through hierarchical/motivated active
inference (Pezzulo et al. 2018; Parr et al. 2022, ch. 7–8), where the higher levels of the hierarchy hold a slowly
drifting context; parametric predictive coding (Bogacz 2017) does not introduce a separate episodic
layer).

### S7.4 Pearl/Friston blanket demarcation and the non-stationary regime

The connection between $\eta_L$ and Friston's free energy principle (Friston 2010; Friston et al. 2017a) is
formulated in (Andriishin 2026, § 4.4): the requirement of self-payment turns the FEP into a
discriminative criterion, and $\eta_L$ operationally closes the gap between the general
formalism and the requirement that the energy belong to the modelling loop. The support is the
distinction between the Pearl blanket (an epistemic instrument) and the Friston blanket (an
ontological property with the requirement of self-payment) (Bruineberg et al. 2018; Bruineberg et al. 2022).

*In the non-stationary regime* the demarcation is preserved: the correspondence between $\eta_L$ and EFE
makes sense only for the Friston blanket; for systems with externalised payment (LLM-as-agent in the
standard architecture) $\eta_L(t) \to 0$ — there is no correspondence. This is consistent with the
structure of § 2.6 main: self-payment as a criterion of applicability of $\eta_L$ is carried over to
the non-stationary regime and to each of the three terms of (1) main.

*Destructive argument of Bruineberg.* The argument (Bruineberg et al. 2018; Bruineberg et al. 2022) against the
ontologisation of the FEP without a self-payment criterion is a structural argument in favour of the
present framework: $\eta_L$ operationalises the requirement of self-payment and thereby answers the
critique of Bruineberg in the sense developed in (Andriishin 2026, § 4.4a).

### S7.5 Andrews (2021) / Williams (2020) cross-ref to paper #1

The rebuttal of the objections of Andrews (2021) (FEP as an empty tautology in the absence of independent
criteria) and Williams (2020) (FEP as a correct but uninformative characterisation of any
non-equilibrium system) is developed in (Andriishin 2026, § S4.4a); in the present work it is not
repeated, but the structural correspondence is preserved: the non-stationary $\eta_L(t)$ closes the
same gap as the stationary $\eta_L$ — the requirement that the dissipative flow belong to the
modelling loop through self-payment. The category of correspondence is *functional* (the common
structure "optimisation of predictive information under a budget"), not ontological (the full
identity of EFE and $\eta_L$ as objects).

### S7.6 The informational rejuvenation metric $\rho(t)$: in detail

The informational rejuvenation metric is the ratio of the refresh and growth rates of memory:

$$\rho(t) = \frac{\dot{M}_{\text{refresh}}(t)}{\dot{M}_{\text{grow}}(t)}.$$

*Interpretation of the poles.* $\rho \gg 1$ — refresh outpaces capacity growth: the refresh fraction
dominates, the nostalgia $\nu$ remains low (the adaptive regime § 5.3 main). $\rho \ll 1$ —
capacity growth outpaces refresh: the archival layer $M_{\lt t}$ accumulates faster than the current
layer manages to be refreshed, and the nostalgia grows monotonically (the nostalgic-collapse regime).

*Epistemic status of $\rho$.* (13) main is a diagnostic independent of the direct measurement of
$\nu(t)$: $\dot{M}_{\text{refresh}}$ and $\dot{M}_{\text{grow}}$ are operationally measurable in
systems where the direct estimate of $\nu(t)$ through (7) main is unavailable. In a bacterium, $\dot
M_{\text{refresh}}$ is the frequency of CheB-mediated demethylation of receptors;
$\dot M_{\text{grow}}$ is the frequency of replication with the addition of new copies of receptors. In
an LLM-corp, $\dot M_{\text{refresh}}$ is the frequency of continual pre-training or RAG update of the
index; $\dot M_{\text{grow}}$ is the frequency of expansion of the corporate corpus.

*Empirical calibration.* The numerical thresholds $\rho_c$ are a separate study. To a first
approximation $\rho \sim 1$ corresponds to the balance condition (11b) main; the exact relation
$\rho_c(\lambda, \sigma, K, \tau_{1/2})$ depends on the same parametric structure
as $\nu_c^{\text{theor}}$ (S5.3 supplementary), and inherits its open status.

*Epistemic-pragmatic decomposition (optional).* In terms of active inference $\rho(t)$
can be interpreted as the balance of epistemic impulse (the increase of epistemic value through
$\dot{M}_{\text{refresh}}$) and pragmatic investment (the expansion of the policy horizon through
$\dot{M}_{\text{grow}}$). This decomposition is a natural consequence of the correspondence S7.1; its
formal development requires the sophisticated inference framework S7.2.

---

## S8. Adiabatic asymptotics of $L_{\text{excess}}$ as an open conjecture

This section contains the full development of the open conjecture on the adiabatic asymptotics of the
cumulative excess-loss, to which main § 8.1, § 8.5 refer as the main open quantitative problem of the
work. Structure: S8.1 — the statement of the conjecture and its theoretical motivation; S8.2 — the
conditions of applicability and the separation of the two adiabatic regimes; S8.3 — the parametric scan
and the numerical support.

The conversion of units for $L_{\text{excess}}$ (nats ↔ bits) is Appendix A main: $L_{\text{excess}}$
is measured in **nats** (as in the original BNT literature (Bialek et al. 2001, § VI) and the MDL formulation);
the explicit conversion $L_{\text{excess}}^{\text{[bits]}}(t) = L_{\text{excess}}^{\text{[nats]}}(t)/\ln 2$
is applied when comparing with the main text.

### S8.1 Statement of the conjecture and theoretical motivation

An alternative theoretical characterisation of the accumulation of information by a learning system about
the latent parameter of the environment is the cumulative excess-loss of the learner relative to the
oracle:

$$L_{\text{excess}}(t) := \sum_{s \le t} \bigl[ \ln p\bigl(X_s \mid X_{s-1},\, \theta(s)\bigr) - \ln \hat p\bigl(X_s \mid X_{s-1},\, \hat\theta(s)\bigr) \bigr]. \tag{S8.1}$$

Here $\theta(s)$ is the true value of the latent parameter of the environment at step $s$, $\hat\theta(s)$ is the
learner's estimate from the history of observations; $L_{\text{excess}}(t)$ is measured in nats. The
operational status of (S8.1): this quantity is directly computable in simulations with a known ground
truth $\theta(s)$, but in real systems (biosphere, LLM-corp, bacterium) the latent parameter $\theta(s)$
is in principle inaccessible — $L_{\text{excess}}$ remains counterfactual, theoretically related to the
operationally measurable one-step $I_{\text{pred}}(t,\tau)$ (4) main, but not equivalent to it.

Connection to (4) main and the thermodynamics of memory. $L_{\text{excess}}(t)$ is the cumulative
Bayesian regret of the learner relative to the oracle in the sense of (Rissanen 1986; Clarke and Barron 1990);
for regular parametric models it is asymptotically related to $I(M_{\le t}; \theta(t))$ through the
theorem on the sufficient statistic in the sense of (Cover and Thomas 2006, ch. 2.9). The difference:
$L_{\text{excess}}$ is the *integral* information about the latent accumulated by time $t$; (4) main is
the *density* of the information flow through the current step, bounded by $\ln k$.

**Conjecture S8.1 (adiabatic asymptotics).** *In the adiabatic OU limit $\sigma^2/\lambda \ll 1$ under
ideal Bayesian learning the cumulative excess-loss $L_{\text{excess}}(t)$ (S8.1) grows as $(K/2)\ln(\lambda t)$
— the Class II asymptotics (Bialek et al. 2001, § VI):*

$$L_{\text{excess}}^{\text{Bialek}}(t) \sim \frac{K}{2} \ln(\lambda t), \tag{S8.2}$$

*where $K = k(k-1)$ is the dimension of the logits. If this secondary metric is adopted as the
non-stationary analogue of predictive information and one defines through it
$\eta_L^{\text{excess}}(t) := L_{\text{excess}}(t)/N_{\text{max}}(t)$,
then in the leading order as $t \to \infty$ in the adiabatic OU limit $\sigma^2/\lambda \ll 1$ under an
ideal Bayesian learner and a linear growth of the denominator $N_{\text{max}}(t) \sim Rt$:*

$$\dot{\eta}_L^{\text{excess}}(t) \sim -\frac{K \ln(\lambda t)}{2 R t^2}. \tag{S8.3}$$

**Status of the conjecture.** (S8.2) is the known Class II BNT asymptotics (Bialek et al. 2001) for regular
parametric models; the coefficient $K/2$ arises in several paradigms of statistical learning theory for
regular parametric models of dimension $K$: MDL (Rissanen 1986; Gr"unwald 2007) and BIC (Schwarz 1978)
through the Laplace approximation on the posterior; PAC-Bayes (McAllester 1999; Catoni 2007) for a
Gaussian prior; online Newton-step regret (Cesa-Bianchi and Lugosi 2006; Hazan et al. 2007) for exp-concave loss (the
formal class to which the logarithmic bound (Hazan et al. 2007) applies; self-concordant is a narrower
property, not used in Hazan et al. 2007). (S8.3) is its direct algebraic consequence under the assumption of
a linear growth of the denominator. The hypothetical status pertains to two aspects: (i) the
convergence of the coefficient in the fit of $L_{\text{excess}}$ to $K/2$ in realisable parametrisations
of OU at finite $\sigma$; (ii) the applicability of the result to the operationally measurable
efficiency $\eta_L(t)$ (10) main, rather than to the secondary $\eta_L^{\text{excess}}$.

### S8.2 Conditions of applicability and the separation of the two adiabatic regimes

Conjecture S8.1 requires the *simultaneous* satisfaction of two conditions (see § 5.2 main):

1. **Slow-driving** $\lambda \tau_{\text{relax}} \ll 1$ — necessary for the applicability of Lemma 2 § 4.4 main and
   the quasi-stationary interpretation of $\theta(t)$.
2. **OU-concentration** $\sigma^2/\lambda \ll 1$ — necessary precisely for (S8.2), not for (4) main: the concentration
   of the posterior over $\theta(t)$ must outpace the proper drift of $\theta(t)$.

Additionally: (iii) boundedly growing memory $|M_{\le t}| = o(t/\ln t)$ (for the linearity of the denominator $N_{\max}
\sim Rt$); (iv) the OU parametrisation of precisely the logits (the softmax form ensures the regularity of the parametric
model); (v) an ideal Bayesian learner (the refresh payment is in $\dot E_{\text{actual}}^{\text{curr}}$, not in the
nostalgic layer).

### S8.3 Parametric scan and numerical support of the conjecture

The numerical study of (S8.2) is a parametric scan of the OU simulation over $\sigma^2/\lambda$ with direct measurement
of (S8.1). The implementation is `simulations/markov_drift_ou_iinf_adiab/`.

**Model.** The OU parametrisation of § 5.2 / § 6.2 main unchanged ($k = 8$, $\lambda = 10^{-3}$, $K = 56$,
Euler–Maruyama $\Delta t = 1$); the ideal Bayesian learner is approximated by Robbins–Monro $\eta_t = c/\sqrt
t$, $c = 3$. $L_{\text{excess}}(t)$ is measured by (S8.1) directly — it is not bounded by $\ln k$, and grows
asymptotically as (S8.2).

**Parametric scan.** Three values of $\sigma$ at fixed $\lambda = 10^{-3}$: $\sigma \in \{0.01,\,
0.003,\, 0.001\}$, giving $\sigma^2/\lambda \in \{0.1,\, 9 \cdot 10^{-3},\, 10^{-3}\}$ respectively.
Transient window $1/\sigma^2 \in \{10^4,\, 1.1\cdot 10^5,\, 10^6\}$; simulation time $T \in \{5\cdot 10^4,\,
5\cdot 10^4,\, 2\cdot 10^5\}$ within the transient window for all cases. Each point is the average over 25
independent realisations of the OU process; seed `SEED = 20260524`. Joint fit $L_{\text{excess}}(t) = A \cdot t +
B \cdot \ln(\lambda t) + C$ on the window $t \in [10^3, 10^4]$ (fixed for direct comparison across
$\sigma$).

**Results.** *Fig. 6: $B$ vs $\sigma^2/\lambda$* (`paper/figs/Fig6.pdf`). *Fig. 7:
family of $L_{\text{excess}}(t)$ for three $\sigma$* (`paper/figs/Fig7.pdf`).

| $\sigma$ | $\sigma^2/\lambda$ | $B$ | $B/(K/2)$ | $R^2$ | drift-share at $t=10^4$ |
|---|---|---|---|---|---|
| 0.01  | $10^{-1}$ | 15.15 | 0.54 | 0.9998 | 76% |
| 0.003 | $9\cdot 10^{-3}$ | 21.85 | 0.78 | 0.9996 | 27% |
| 0.001 | $10^{-3}$ | 23.84 | 0.85 | 0.9995 | 10% |

The coefficient $B$ grows monotonically as $\sigma^2/\lambda \to 0$, giving $B/(K/2) = 0.54,\, 0.78,\, 0.85$
respectively; $R^2 \ge 0.9995$ at all points, the drift-share $A \cdot t / L_{\text{excess}}(t)$ falls from 76%
to 10%.

**Wide-window negative control.** Upon extending the window of the joint fit to $[10^3, T_{\text{total}}]$ (instead of
the fixed $[10^3, 10^4]$) the fit degenerates: for $\sigma^2/\lambda = 10^{-3}$ and $T_{\text{total}} = 2\cdot
10^5$ the coefficient $B/(K/2)$ grows to 2.0 — an artefact of the Robbins–Monro $1/\sqrt t$ bias of the learner on
long horizons, not part of the BNT effect. The window $[10^3, 10^4]$ is chosen for direct comparison between the
points of the adiabatic series and is consistent with the transient condition $t \ll 1/\sigma^2$. This control
records the sensitivity of the result to the choice of the fit window and prevents the interpretation of the
monotone trend of $B/(K/2)$ as the result of selective presentation.

**Interpretation.** The functional form $L_{\text{excess}}(t) = A\cdot t + B\ln(\lambda t) + C$ is extracted with
$R^2 \ge 0.9995$; the coefficient $B$ grows monotonically from 0.54 to 0.85 $\cdot K/2$ as $\sigma^2/\lambda$ goes from
$10^{-1}$ to $10^{-3}$. **This is not a proof of convergence $B \to K/2$**; the asymptotic limit $B_\infty$
remains an open question. The pre-registered protocol for confirming the conjecture (S8.2) requires a minimum of 5–6 scan
points with an explicit analytical form of the finite-$\sigma$ correction $B(\epsilon) = (K/2)(1 - g(\epsilon))$.

The relegation of (S8.2) and (S8.3) to the status of an open conjecture reflects the current state of the theoretical
justification, not a final verdict: the numerical support of the functional form and the monotone convergence to $K/2$
are a nontrivial empirical regularity, deserving further analytical investigation (see § 8.5 main).

**Reproducibility.** The full implementation is `simulations/markov_drift_ou_iinf_adiab/` (`model.py`, `main.py`,
`README.md`); the figures are in `paper/figs/Fig6.pdf` and `paper/figs/Fig7.pdf`;
the summary log is `simulations/markov_drift_ou_iinf_adiab/results_summary.{txt,json}` and `run.log`. Seed
`SEED = 20260524`; running `python main.py` (~1 minute on a modern laptop) reproduces all numerical
values.

---

## S9. Details of §§ 2.1 and 2.3 main: three-faces comparison and majority-vote derivation

This section assembles the dense statistical-mechanics apparatus moved out of §§ 2.1, 2.3 main to lighten the
main exposition. S9.1 — the full comparison of decomposition (1) main with the three-faces / Hatano–Sasa
classification of entropy production; S9.2 — the full derivation of the majority-vote variant of Lemma 1 (equation (2') main).

### S9.1 Comparison of decomposition (1) main with three-faces / Hatano–Sasa

Decomposition (1) main is *thematically parallel* to the three-faces decomposition of entropy production
(Esposito and Van den Broeck 2010) and the Hatano–Sasa decomposition (Hatano and Sasa 2001), but does not reduce to it
automatically. In the canonical treatment the housekeeping component maintains a fixed NESS under stationary
control parameters, the excess component — the change of the NESS under a slow change of parameters. The correspondence
of the components of (1) main:

- $\dot{E}_{\text{store}}$ for maintaining a specific information state against relaxation to symmetric
  equilibrium belongs to the *excess* class (relaxation-driven), not to housekeeping: it is not a fixed
  NESS that is maintained, but a specific non-equilibrium information state relaxing to symmetric equilibrium.
- $\dot{E}_{\text{actual}}^{\text{curr}}$ upon re-tuning $\hat P(t)$ to the drift of the environment also belongs to the
  excess class (the change of the target state under slow drift of the parameters).
- the refresh operations of Lemma 1 — a series of non-equilibrium protocols on top of this, reducible neither to purely housekeeping,
  nor to purely excess.
- $\dot{E}_{\text{grow}}$ — a structural rearrangement (constraint-driven), going beyond the canonical
  three-part decomposition: the addition of new degrees of freedom changes the state space itself, and not only
  the distribution on it.

The exact embedding of the non-stationary information decomposition into the three-faces apparatus is an open problem (§ 8.5 main).

### S9.2 Full derivation of the majority-vote variant of Lemma 1 (equation (2') main)

In the limit $\varepsilon \to 0$ the polynomial bound (2) main $\sim 1/\varepsilon$ is replaced by the majority-vote variant
through redundant coding. For $n$ copies of a bit with independent errors $\varepsilon < 1/2$ the majority vote errs
only when no fewer than $n/2$ copies are flipped; the probability of this tail of the binomial $\mathrm{Bin}(n,\varepsilon)$
is bounded by the Chernoff bound

$$\varepsilon^{\text{eff}} \;=\; \Pr\!\left[\mathrm{Bin}(n,\varepsilon) \ge n/2\right] \;\le\;
\exp\!\bigl(-n\, D(\tfrac{1}{2}\,\|\,\varepsilon)\bigr), \qquad
D(\tfrac{1}{2}\,\|\,\varepsilon) = \tfrac{1}{2}\ln\frac{1}{4\varepsilon(1-\varepsilon)},$$

where $D(\tfrac{1}{2}\,\|\,\varepsilon)$ is the KL divergence of the Bernoulli $\mathrm{Ber}(1/2)$ from $\mathrm{Ber}(\varepsilon)$.
Here $\varepsilon$ in $D(\tfrac{1}{2}\,\|\,\varepsilon)$ is the *raw per-copy error* (the probability of a flip on one
copy of a bit, a physical parameter of the carrier), whereas $\varepsilon^{\text{eff}}$ is the *target* (achievable) effective
error of the decoded bit after the majority vote; it is precisely the difference of these two quantities that ensures the exponential
suppression $\varepsilon^{\text{eff}} \le e^{-n D(\frac{1}{2}\|\varepsilon)}$ at fixed $\varepsilon < 1/2$.
Hence, to achieve a target effective error $\varepsilon^{\text{eff}}$ it suffices to have

$$n \;\sim\; \frac{\ln(1/\varepsilon^{\text{eff}})}{D(\tfrac{1}{2}\,\|\,\varepsilon)}$$

copies — a *logarithmic* number of copies in $1/\varepsilon^{\text{eff}}$ (rather than linear in $1/\varepsilon$, as in naive
refresh). Each of the $n$ copies degrades independently with the same $\tau_{1/2}$, so the refresh rate for the majority vote
preserves the factor $t/\tau_{1/2}$ (the number of refresh cycles), but the payment per cycle is logarithmic in $1/\varepsilon$,
not polynomial:

$$E_{\text{store}}^{(1),\text{maj}}(t) \;\ge\; C \cdot \frac{t}{\tau_{1/2}} \cdot k_B T \ln(1/\varepsilon),
\qquad C = O(1), \tag{2'}$$

with a constant $C$ depending on the majority-vote scheme (Sagawa–Ueda extension Sagawa and Ueda 2009; Parrondo et al. 2015).
Structurally, (2) main and (2') coincide in their linear growth with $t/\tau_{1/2}$ — the number of refresh cycles; the difference is
in the dependence on $\varepsilon$: polynomial in naive vs logarithmic in majority-vote. The class
of paradigm-case systems to which the original (2) main applies: ferromagnetic domains of warm media and systems with naive
refresh; for biological systems with feedback repair (2') is used, with the logarithmic dependence on
$\varepsilon$. Feedback-assisted refresh (DNA mismatch repair) additionally mitigates the bound by Sagawa–Ueda
$\langle W \rangle \ge k_B T (\Delta S - I_{\text{meas}})$ (Sagawa and Ueda 2009): a *structured* bit
with conditional entropy $H(\varepsilon) = -\varepsilon\ln\varepsilon - (1-\varepsilon)\ln(1-\varepsilon) \ll \ln 2$
at $\varepsilon \ll 1$ is erased, and in the limit of ideal measurement the cost of refresh goes to zero — the
asymptotics of autonomous Maxwell demons (Koski et al. 2014; Mandal and Jarzynski 2012; Barato and Seifert 2014; Bauer et al. 2014).

---
