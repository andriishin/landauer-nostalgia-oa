<!--
  English version of the article (submission language for Theory in Biosciences).
  Translated from main.ru.md (primary source).
  Translation glossary — docs/translation-notes.md (RU→EN map for central terms;
  "informational nostalgia" as the mandatory EN form of the central term).
-->

# The Non-Stationary Landauer Efficiency: Memory Growth, Informational Nostalgia, and the Thermodynamics of Learning Systems

**Authors:** Alexander Andriishin ¹\*

¹ Independent researcher, ORCID: 0009-0000-4739-4017 (https://orcid.org/0009-0000-4739-4017)
\* Correspondence: alexander@andriishin.com

## Abstract

Memory is not a record but an active thermodynamic process that continuously pays for itself. A bacterium holds a ligand concentration through receptor methylation; a mammalian embryo erases its epigenetic past through reprogramming; a continual-learning network spends energy updating its weights; a biosphere keeps in its genomes a chronicle of a vanished environment. In each, memory is held against an environment that drifts faster than it is updated, and a growing fraction of the stored bits no longer predicts anything. We call this fraction *informational nostalgia* (a technical term, not a psychological one) and extend the Landauer efficiency to the non-stationary regime under one condition: every bit must be paid for by the system's own free energy (the self-payment postulate). Two results are proved. *Lemma 1*: holding a bit against thermodynamic erosion requires strictly positive power, with an explicit lower bound. *Lemma 2*: for any finite refresh rate under slow Ornstein–Uhlenbeck drift, the nostalgic fraction asymptotically stays above an explicit positive constant — staleness is inevitable. Two consequences follow: efficiency vanishes above a critical nostalgia, and an optimal rate of complete reset exists, of an order matching epigenetic reprogramming in mammalian embryogenesis and seasonal diapause in *Daphnia*. The framework binds bet-hedging, active inference, and catastrophic forgetting within one thermodynamic budget; a pre-registered *E. coli* chemotaxis prediction makes it falsifiable on a standard protocol. The stationary case is developed in a companion work on the Landauer efficiency of self-modeling and recovered here as a strict limiting case.

## Keywords

informational nostalgia; Landauer's principle; non-stationary thermodynamics; biological memory;
active inference; bet-hedging

---

## 1. Introduction

*E. coli chemotaxis as a motivating example.* The bacterium *Escherichia coli* swims up an aspartate gradient, holding the memory of the ligand concentration a minute ago through the covalent methylation of four residues of the Tar/Tsr receptors (Sourjik and Berg 2002). The methyltransferase CheR adds methyl marks, the methylesterase CheB removes them; the equilibrium between them is precisely the "memory" of the recent gradient, translated into a tumbling bias. The characteristic refresh time of these marks is on the order of a minute, the response time tens of seconds. When the aspartate source moves abruptly, the methyl marks update more slowly than the gradient has shifted: part of the receptor configuration encodes an environment that no longer exists. This *stale fraction of memory* keeps costing the bacterium bits of ATP (maintenance of methyltransferase activity, repair of degraded marks) but does not help predict the local gradient now. That sensory adaptation in this system is coupled to an energetic cost is an established fact: the energy–speed–accuracy bound for chemotaxis adaptation (Lan et al. 2012) and the general energetics of cellular computation (Mehta and Schwab 2012) showed that the accuracy and rate of adaptation are paid for by dissipation. Adaptation catches up — but not instantaneously; the lag between an environmental shift and a memory shift is physically costly. From this follows a concrete falsifiable outcome: a bacterium with reduced methyltransferase activity ($cheR^{\downarrow}$, a slowed mark-refresh rate) under periodic switching of attractant concentration should recover adaptation *more slowly* than the wild type — a measurable slowdown of the adaptation rate on a standard laboratory protocol (§ 8.6, Prediction 2).

*The same structure recurs across scales.* Biospheres accumulate genomic content through evolutionary transitions (Maynard Smith and Szathm'ary 1995; Landenmark et al. 2015); abiotic cycles turn over on aeonic windows, leaving in genomes a chronicle of past regimes. Populations in fluctuating environments maintain phenotypic diversity as an informational bet on possible future states (Kussell and Leibler 2005) — bet-hedging describes the same asymmetry between the refresh rate of the phenotype distribution and the rate of the environmental shift. In mammals, early embryogenesis features the almost complete erasure of DNA methylation followed by its re-establishment (Reik 2007) — a massive reset of stale epigenetic memory. Neural networks under continual learning exhibit *catastrophic forgetting* (McCloskey and Cohen 1989; Kirkpatrick et al. 2017); large language models with weights frozen after training accumulate ever more parameters describing the vocabulary and facts current at the moment of cutoff but drifting with the stream of new data (Lazaridou et al. 2021). All these systems share a common structural feature: the memory of a past environmental state grows, the environment drifts on the same timescales, and the optimal strategy of memory retention balances the informational gain against the thermodynamic cost. A quantitative measure of this balance is an open question at the intersection of stochastic thermodynamics, population biology, and learning theory.

*The measure of the balance.* To measure this balance one needs a single dimensionless quantity that sets the informational gain of memory against the thermodynamic cost of its maintenance. Such a measure is the *Landauer efficiency* $\eta_L$ — intuitively, the fraction of the thermodynamic payment over the entire history of memory that is actually convertible into prediction of the environment. Formally,

$$\eta_L(t) = \frac{I_{\text{pred}}(t,\tau)}{N_{\text{max}}(t)},$$

where the numerator $I_{\text{pred}}(t,\tau)$ is the predictive information (Tishby et al. 1999; Bialek et al. 2001), the bits held by the system's memory about the prediction horizon $\tau$ of its own environment, and the denominator $N_{\text{max}}(t)$ is the cumulative Landauer budget: all the free energy spent by the system up to time $t$, expressed in bits through the divisor $k_B T \ln 2$. In the stationary formulation (Andriishin 2026) the budget is the operational payment for a single window; in the non-stationary regime it accumulates with history and decomposes into several qualitatively distinct flows (detailed definition — § 2.1). The efficiency is bounded above, $\eta_L \le 1$, only when the numerator and denominator belong to the same physical system: the *self-payment postulate* — the dissipative flow must be paid for by the system's own metabolism rather than externalised, since the system cannot convert into prediction more bits than it has itself paid for with its own free energy; energy borrowed from outside would lift the ceiling. The operational test is counterfactual ablation (cutting off the metabolic supply and checking the degradation of memory; § 2.6, § 4.1). The stationary case $\dot E_{\text{store}} = \dot E_{\text{grow}} = 0$ is developed in (Andriishin 2026); the present work extends the framework without refuting it.

*Content of the work.* Two structural results are proved, with an intuitive interpretation:

- **Lemma 1** (§ 2.3) — every held bit of memory costs nonzero power (for feedback memory — in form (2')); concretely,
  $\dot E_{\text{store}}^{(1)} \ge (\ln 2)/(2\varepsilon\tau_{1/2}) \cdot k_B T \ln 2$ for the
  naive-refresh scheme at error threshold $\varepsilon$ and accuracy half-life
  $\tau_{1/2}$ (ferromagnetic domains, neural-network weights without a correlated reference degrade);
  for feedback memory with a correlated reference (methylation, epigenetic marks —
  CheR/CheB repair) — the form (2') of feedback-assisted refresh, logarithmic in $\varepsilon$.
- **Lemma 2** (§ 4.4) — under slow OU drift of the environment and finite $\dot M_{\text{refresh}}$
  the nostalgic fraction $\nu^{\text{theor}}(t)$ asymptotically does not fall below $1 - c$,
  where $c = \dot M_{\text{refresh}} \tau_E / \lvert M_{\le t}\rvert$ is an explicit computable
  constant from the refresh rate and the characteristic drift time.
- **The critical value $\nu_c$** (§ 5.3) — there exists a nostalgia threshold above which
  $\eta_L(t)$ vanishes; for systems with a fixed balance ratio $\eta_M$.
- **The optimal reset $\tau_{\text{reset}}^*$** (§ 5.3) — there exists an optimal rate of
  complete memory reset (biological example: epigenetic reprogramming in early mammalian
  embryogenesis (Reik 2007); *Daphnia* diapause on seasonal windows).

*Distinction from adjacent frameworks.* The parallels with bet-hedging, catastrophic forgetting, and concept drift are structural, not equivalences: what is new here is the common thermodynamic scale on which the numerator and denominator of the efficiency are paid for by one physical system, and the lower bound on inevitable nostalgia (Lemma 2). A detailed comparison — § 8.4.

*Structure of the work.* § 2 — the formalism of the budget and self-payment, glossary of key notation (§ 2.0); § 3 — predictive information on growing memory; § 4 — nostalgia and Lemma 2; § 5 — $\eta_L(t)$, three regimes, the regime threshold $\nu_c$, the optimal $\tau_{\text{reset}}^*$; § 6 — numerical illustration; § 7 — connection with active inference; § 8 — summary, applied consequences, open problems.

## 2. The Non-Stationary Landauer Budget $N_{\text{max}}(t)$

### 2.0 Glossary of key notation and abbreviations

An extended notation table — Supplementary § S3; below is the minimal set sufficient for reading § 2–5 without consulting the supplementary.

*Abbreviations.*

- **OU** — Ornstein–Uhlenbeck stochastic process: a Gaussian process with exponential
  autocorrelation; in this work — the canonical class of "slow environmental drift" (§ 5.2).
- **DPI** — Data Processing Inequality: $I(X; Y) \ge I(X; f(Y))$ for any processing $f$
  (Cover and Thomas 2006, § 2.8).
- **KSG** — Kraskov–Stögbauer–Grassberger MI estimator: a $k$-NN estimate of mutual information from
  samples (Kraskov et al. 2004).
- **MI** — mutual information.
- **EFE** — expected free energy (Friston et al. 2017a); see § 7.
- **FIFO** — first-in-first-out memory refresh scheme; § 4.4.
- **PSP** — piecewise-stationary process (a piecewise-stationary process with Poisson regime switches);
  numerical illustration of drift § 6.1.
- **$\nu^{\text{theor}}$ vs $\nu^{\text{op}}$** — *theoretical* nostalgia (via a
  counterfactual oracle with access to the true $\theta(t)$) and *operational* nostalgia (via a
  measurable proxy $C_\mu^{\text{emp}}$); related by the DPI inequality
  $\nu^{\text{op}} \le \nu^{\text{theor}}$ (§ 4.1).

*Summary table of timescales* (six quantities):

| Symbol | Meaning | Biological example |
|---|---|---|
| $\tau_E$ | characteristic environmental drift time (OU correlation time, $\tau_E = 1/\lambda$) | minutes for ligand gradients; years for language vocabulary; aeons for biospheric cycles |
| $\tau_{1/2}$ | half-life of the *excess* accuracy of reading a bit ($p - 1/2$, not $p$ itself; a symmetric channel relaxes to $p \to 1/2$) | seconds–minutes for receptor methyl marks; years for LLM weights under drift during inference |
| $\tau_{\text{relax}}$ | relaxation of the chain to the stationary distribution | $\tau_{\text{relax}} = 1/(2\gamma)$ for the symmetric two-state channel of Lemma 1 |
| $\tau_d$ | characteristic response time of the system (window of counterfactual ablation) | $\sim 30$ min for a bacterium; $\sim 10$ years for the biosphere (response time, not aeon); a quarter for a corporation |
| $\tau_{\text{reset}}^*$ | optimal interval of periodic memory reset | $\sim$generation time for epigenetic reprogramming (Reik 2007); $\sim$season for *Daphnia* diapause |
| $\tau_w$ | sliding window of the empirical estimate $\eta_L^{(\tau_w)}(t)$ | the experimenter's choice; a parameter of the analysis, not of the system |

Detailed definitions with dimensions, value ranges, and defining formulas — Supplementary § S3.1–S3.9.

### 2.1 Sources of dissipation in the non-stationary regime

In the stationary formulation (Andriishin 2026, § 2.1) the denominator of $\eta_L$ is defined as $N_{\text{max}} = E_{\text{actual}} / (k_B T \ln 2)$, where $E_{\text{actual}}$ is the total exergy that can physically pay for Landauer erasures over the characteristic window. This definition implicitly relies on stationarity: exergy and predictive information are collected over one and the same window $\tau$ and normalised on it; the temporal profile within the window is immaterial.

In the non-stationary regime this assumption is violated. The system's total payment for maintaining its own model from $t = 0$ to the present moment $t$ is composed of three qualitatively distinct components. The first is the current operational dissipation $E_{\text{actual}}^{\text{curr}}(t)$, paying for the incremental erasure and updating of the model on the interval $[t - \Delta t,\, t]$; structurally it coincides with the stationary $E_{\text{actual}}$ but has meaning only as instantaneous power. The second is the cumulative cost of holding old information:

$$E_{\text{store}}(t) = \int_0^t \dot{E}_{\text{store}}(t')\,dt',$$

where $\dot{E}_{\text{store}}(t')$ is the power required to maintain the already accumulated bits of memory $M_{\lt t'}$ under conditions of thermodynamic erosion. Without a supply of free energy any physical memory degrades — nucleic acid hydrolyses, epigenetic marks are lost upon replication, neural-network weights are blurred by numerical noise during long-term storage on flash media. Holding information against this erosion has a physical cost. The third is the cost of memory expansion:

$$E_{\text{grow}}(t) = \int_0^t \dot{E}_{\text{grow}}(t')\,dt',$$

where $\dot{E}_{\text{grow}}(t')$ is the power that physically pays for creating new memory capacity. In a biological system this is the synthesis of biopolymers for DNA replication and the formation of new synaptic contacts; in a learning artificial system — weight-write operations; in a city — the building of archives, infrastructure, and educational institutions.

The full cumulative Landauer budget takes the form

$$N_{\text{max}}(t) = \frac{1}{k_B T \ln 2} \left( \int_0^t \dot{E}_{\text{actual}}^{\text{curr}}(t')\,dt' + E_{\text{store}}(t) + E_{\text{grow}}(t) \right). \tag{1}$$

The denominator of $\eta_L(t)$ is precisely $N_{\text{max}}(t)$ by (1) — a cumulative quantity, not a sliding average: in the non-stationary case the cumulative nature of $N_{\text{max}}(t)$ retains information about the entire history of the payment and ensures a correct comparison with the predictive information accumulated up to $t$. The decomposition (1) is *thematically parallel* to the three-faces decomposition of entropy production (Esposito and Van den Broeck 2010) and the Hatano–Sasa decomposition (Hatano and Sasa 2001), but does not reduce to it automatically: $\dot{E}_{\text{store}}$ and $\dot{E}_{\text{actual}}^{\text{curr}}$ belong to the *excess* class (relaxation-driven), the refresh operations of Lemma 1 are a series of non-equilibrium protocols layered on top of it, and $\dot{E}_{\text{grow}}$ (constraint-driven) falls outside the canonical three-part decomposition; the exact embedding of the non-stationary information decomposition into the three-faces apparatus is an open problem (§ 8.5). A detailed comparison with the housekeeping/excess classification — Supplementary § S9.1. The alternative through a sliding window (§ 5.1) — the instantaneous efficiency $\eta_L^{(\tau)}(t)$, which has independent meaning but does not close the problem (Andriishin 2026, § 6).

### 2.2 Consistency with the stationary case

*The stationary-accounting postulate (self-contained).* The postulate on which the conditional status of the bound $\eta_L \le 1$ rests asserts the following: in the stationary regime (zero growth and refresh of memory) the system's entire information accounting passes through a single flow — the payment for storage $\dot{E}_{\text{actual}}^{\text{curr}}$ — and the bound $\eta_L \le 1$ holds under the same two premises as in (Andriishin 2026): self-payment by the system for its own model, and the reducibility of the entire payment to this single flow. Below, this consistency is traced with navigation to specific steps of paper #1.

In the stationary regime with $\dot{M} = 0$ we have $\dot{E}_{\text{grow}} = 0$; under optimal forgetting $\dot{E}_{\text{store}}$ reduces to a refresh of the current content and enters $\dot{E}_{\text{actual}}^{\text{curr}}$. (1) reduces to $N_{\text{max}}(\tau) = E_{\text{actual}}^{\text{curr}} \tau / (k_B T \ln 2)$ (Andriishin 2026, (1) § 2.1); the stationary lemma $E_{\text{actual}} \ge I_{\text{pred}} \cdot k_B T \ln 2$ is a special case at $\dot{E}_{\text{store}} = \dot{E}_{\text{grow}} = 0$ **and given the truth of the stationary-accounting postulate of paper #1** (Andriishin2026, § 2.1 step (iv) of the Lemma + condition (c) of the Proposition). This postulate asserts that in the stationary regime, with zero memory growth and refresh, the system's entire information accounting reduces to the storage balance (the entire payment passes through $\dot{E}_{\text{actual}}^{\text{curr}}$); it is precisely this on which the conditional status of the bound $\eta_L \le 1$ depends. Reducibility to the stationary case is ensured by two conditions on the budget flows ($\dot{E}_{\text{store}} = \dot{E}_{\text{grow}} = 0$), does not require the limiting transition $\tau_E \to \infty$ in the characteristic environmental drift time; moreover, the one-step predictive information $I_{\text{pred}}^{\text{stat}} = I(X_t; X_{t+\tau} \mid \theta_0)$ (§ 3.1) is a finite constant of the ergodic chain. Under the additional conditions of first-order Markovianity with $M_t = X_t$, $\hat P(t) \to P_0$, and the choice $\tau = \tau_d$ of paper #1 (by the convention of § 2.1 Andriishin 2026), the one-step $I_{\text{pred}}$ coincides with the $I_{\text{pred}}^{\text{stat}}$ of paper #1, and the Lemma on stationary erasure (Andriishin 2026, § 2.1) is **consistent with the stationary lemma as a special case of (1) under the stated conditions and within the domain fixed by the Proposition of paper #1**; the "conditional" status — $\eta_L \le 1$ as a conditional consequence of Landauer's principle given the truth of both supporting premises (self-payment + the stationary-accounting postulate) — is inherited by the non-stationary extension: (1) and the reducibility to the stationary case depend on the same postulate (c); together they constitute a conditional extension of the stationary result, not a theorem proved for an arbitrary refresh channel.

*Difference of operationalisations.* The core definition (4) — the one-step conditional MI between environmental states with horizon $\tau$ — is structurally distinct from the original $I_{\text{pred}}^{\text{stat}}(X, \tau) := I(M_t; X_E^{[t,t+\tau]})$ of paper #1 (Andriishin 2026, (1a) § 2.1), where the predictive information is defined as the MI between the memory configuration and the environmental trajectory over a window. For a first-order Markov chain with $M_t = X_t$ and sufficient $\hat P(t) = P_0$, both quantities yield one constant at $\tau = 1$; in the general case the one-step (4) corresponds to the entropy rate $h_\mu$, and the interval-based one to the excess entropy $C_\mu$ (Crutchfield–Feldman). The numerical estimates (Andriishin 2026, § 3.5) for specific systems carry over as order-of-magnitude references; exact equality is attained only in the stated special case. The $\eta_L^{\text{stat}}$ of the present work and the quantities of (Andriishin 2026, § 2.1) structurally correspond to one another with this caveat.

### 2.3 Lemma on the minimal cost of storage

The central statement of the section relates $\dot{E}_{\text{store}}$ to the Landauer bound through cyclic rewriting.

**Lemma 1 (minimal storage cost, naive refresh).** *Consider a symmetric two-state channel $\{0, 1\}$ with a symmetric transition-rate matrix $\gamma$, stationary distribution $p^* = 1/2$, and master equation $\dot{p} = -2\gamma(p - 1/2)$. Let a bit of information be held in this channel on the interval $[0, t]$ with transition rate $\gamma$ (relaxation time to equilibrium $\tau_{\text{relax}} = 1/(2\gamma)$; the half-life of the* excess *reading accuracy $p - 1/2$ — the channel is symmetric and relaxes to $p \to 1/2$, so it is the deviation from $1/2$ that half-decays, not $p$ itself — $\tau_{1/2} = (\ln 2)/(2\gamma)$). Let it be required to hold the probability of correctly reading the bit no lower than $1 - \varepsilon$ for a fixed threshold $\varepsilon \in (0, 1/2)$, with the restoration procedure being a naive refresh: periodic resetting of the noisy register followed by writing the reference value without using measurement feedback. Since naive refresh resets an* unstructured *bit ($H = \ln 2$) without a correlated reference, the Sagawa–Ueda mitigation (Sagawa and Ueda 2009) is inapplicable to it, and each reset is logically irreversible ($\ge k_B T \ln 2$ by Landauer 1961) — this is a premise of the estimate, not a consequence. Then the minimal dissipation for holding this bit is bounded below by the quantity*

$$\dot E_{\text{store}}^{(1)} \;\ge\; \frac{\ln 2}{2\varepsilon \cdot \tau_{1/2}} \cdot k_B T \ln 2, \qquad
E_{\text{store}}^{(1)}(t) \;\ge\; \frac{t \cdot \ln 2}{2\varepsilon \cdot \tau_{1/2}} \cdot k_B T \ln 2. \tag{2}$$

*Dimensions.* The left-hand side of (2, first) is power (W): $\dot E_{\text{store}}^{(1)}$ is the instantaneous payment for holding a single bit; the right-hand side is $[k_B T] \cdot [1/\tau_{1/2}] = $ J/s. The integral form (2, second) is obtained by multiplying by $t$: $E_{\text{store}}^{(1)}(t)$ in J. Both forms are used below according to context: Eqs. (1) and (8) are flows $\dot E$; the cumulative estimates of § 6 are $E(t)$.

*Relation to the energy–accuracy bound of Lan et al. (2012).* Lemma 1 does not reintroduce energetics into adaptation: the energy–speed–accuracy bound (Lan et al. 2012; Mehta and Schwab 2012) already relates the accuracy of a *single* adaptation to dissipation. The contribution of Lemma 1 is the non-stationary extension of this balance to the *holding* of an ageing bit over time: not the price of reaching the adapted state, but the minimal flow for maintaining the readability of memory against relaxation, accumulating as $E_{\text{store}}^{(1)}(t) \propto t$.

**Proof.** A symmetric two-state Markov channel with transition rate $\gamma$ gives the equation $\dot{p}(t) = -2\gamma (p - 1/2)$, whence $p(t) - 1/2 = (1/2) e^{-2\gamma t}$ under the initial condition $p(0) = 1$. The condition $p(\Delta t) \ge 1 - \varepsilon$ over the inter-refresh interval gives $e^{-2\gamma \Delta t} \ge 1 - 2\varepsilon$, whence $\Delta t \le -\ln(1 - 2\varepsilon)/(2\gamma) \approx \varepsilon/\gamma$ for small $\varepsilon$. The minimal refresh frequency: $\nu_{\text{refresh}} \ge \gamma/\varepsilon = (\ln 2)/(2\varepsilon \cdot \tau_{1/2})$ — which is precisely the right-hand side of (2, first) after multiplying by $k_B T \ln 2$. Over the interval $[0,t]$ at least $t \cdot (\ln 2)/(2\varepsilon \cdot \tau_{1/2})$ refresh operations are performed; each operation includes at least one logically irreversible reset of the noisy register (without a correlated reference) and, by Landauer's theorem (Landauer 1961), dissipates $\ge k_B T \ln 2$ — logical irreversibility and the inapplicability of the Sagawa–Ueda mitigation are fixed by the premise of the lemma. By construction, naive refresh *does not read* the current register state (open-loop, scheduled rewriting), and therefore cannot use the correlation "current register $\leftrightarrow$ reference"; $H = \ln 2$ refers to the content being erased at the moment of reset, regardless of whether the bit has had time to degrade. Integrating $\dot E_{\text{store}}^{(1)}$ over $[0,t]$ gives (2, second). $\square$

*Dependence on $\varepsilon$ and feedback-assisted refresh.* The polynomial factor $1/\varepsilon$ in (2) is a consequence of naive refresh without amortisation; the logarithmic dependence $\sim k_B T \ln(1/\varepsilon)$ is recovered by a block majority-vote strategy ($n \sim \ln(1/\varepsilon^{\text{eff}})/D(\tfrac{1}{2}\|\varepsilon)$ copies (Sagawa and Ueda 2009; Parrondo et al. 2015), requires $E_{\text{grow}}$). For feedback-assisted refresh (DNA mismatch repair) the bound is mitigated, à la Sagawa–Ueda, towards the asymptotics of autonomous Maxwell demons (Koski et al. 2014; Mandal and Jarzynski 2012; Barato and Seifert 2014; Bauer et al. 2014) (full derivation — Supplementary § S9.2). (2) applies to non-feedback memory without a correlated reference (ferromagnetic domains; neural-network weights under the non-feedback treatment); for feedback memory with a correlated reference (methylation/epigenetics — CheR/CheB) — the form (2'), an upper estimate.

*Behaviour in the limit $\varepsilon \to 0$.* In the physically interesting limit $\varepsilon \to 0$ the bound (2) is polynomial $\sim 1/\varepsilon$ and useless for systems with small $\varepsilon$ (biological memory with $\varepsilon \sim 10^{-3}$). In this regime (2) is replaced by a *majority-vote variant* through redundant coding with $n \sim \ln(1/\varepsilon^{\text{eff}})/D(\tfrac{1}{2}\|\varepsilon)$ copies, giving a logarithmic payment per cycle:

$$E_{\text{store}}^{(1),\text{maj}}(t) \;\ge\; C \cdot \frac{t}{\tau_{1/2}} \cdot k_B T \ln(1/\varepsilon),
\qquad C = O(1) \tag{2'}$$

(the Sagawa–Ueda extension (Sagawa and Ueda 2009; Parrondo et al. 2015); full derivation through the Chernoff bound $\varepsilon^{\text{eff}} \le \exp(-n\,D(\tfrac{1}{2}\|\varepsilon))$, whence $n \sim \ln(1/\varepsilon^{\text{eff}})/D(\tfrac{1}{2}\|\varepsilon)$ — Supplementary § S9.2). Structurally (2) and (2') coincide in their linear growth with $t/\tau_{1/2}$; the difference is in the dependence on $\varepsilon$ (polynomial in naive vs logarithmic in majority-vote). The original (2) applies to naive-refresh memory (ferromagnetic domains); for biological systems with feedback repair, (2') is used.

Lemma 1 *complements* the stationary lemma (Andriishin 2026, § 2.1), it does not generalise it: it bounds a new component $\dot{E}_{\text{store}}$ not defined in the stationary regime, under the same domain conditions (including the stationary-accounting postulate of the Proposition § 2.1 of paper #1). Consistency with § 2.2 — through $\dot{E}_{\text{store}} = \dot{E}_{\text{grow}} = 0$, not through $\tau_{1/2} \to \infty$. In real systems $\tau_{1/2}$ is finite (DNA replication, drift of magnetic domains, numerical error of weights), and (2) gives a nonzero price. (2) is a particular application of (Sagawa and Ueda 2009) to the problem of maintaining a bit; the general apparatus is the stochastic thermodynamics of computation (Wolpert 2019; Hartich et al. 2014).

### 2.4 The cost of memory expansion

$E_{\text{grow}}(t)$ pays for creating new capacity — moving a physical degree of freedom from an equilibrium state into a controlled one. This is the inverse of erasure: the work is bounded below by

$$E_{\text{grow}}^{(1)} \ge k_B T \ln 2 \tag{3}$$

per bit of added capacity (the standard isothermal argument — the Szilard cycle; the inverse of the Landauer erasure bound). In real systems the actual cost exceeds (3) by many orders of magnitude: the synthesis of one nucleotide during DNA replication $\sim 2$ molecules of ATP $\approx 10^{-19}$ J — five to six orders of magnitude larger than $k_B T \ln 2 \approx 3 \cdot 10^{-21}$ J/bit at $T = 310$ K. This is consistent with the observation (Andriishin 2026, § 2.1) that $E_{\text{comp}}/E_{\text{actual}} \ll 1$.

### 2.5 Interpretation and limits of applicability

Formula (1) fixes a structural result: accumulated history has a physical price on a par with operational dissipation. For slowly growing memory the stationary approximation is workable; for rapidly growing memory the correction dominates. Limits of applicability: (i) separability of the three flows; (ii) the irreversible computational regime (Bennett 1973; Bennett 1982) — for Bennett-reversible computation $N_{\text{max}}$ loses meaning (Andriishin 2026, § 2.1).

### 2.6 Self-payment in the non-stationary regime

*Why biology should measure self-payment.* Before formulating the criterion, let us explain why it is substantive and not ad hoc. Memory is cheap to write off as "information in the substrate" — but the information in the DNA of a museum specimen is also measurable, even though no loop maintains it any longer. The difference between actively maintained memory and an archival record is precisely that the carrier system *itself pays* for holding its model of the environment: it expends its own free energy against thermodynamic erosion. Therefore the quantity $-dF_X/dt$, measured over the response window $\tau_d$, is not an arbitrary indicator but a direct counter of whether the system maintains its memory actively. The self-payment criterion makes measurable exactly this boundary: $\eta_L$ is meaningful where the numerator and denominator are paid for by one self-payment loop, and loses meaning where the payment is borne by an outside observer.

*The self-payment postulate* (Andriishin 2026, § 2.7): to apply the ratio $\eta_L$ to a physical system, it is required that the numerator ($I_{\text{pred}}$) and the denominator ($N_{\text{max}}$) belong to *one* self-payment loop — the dissipative flow must be paid for by the system's own metabolism rather than externalised onto an outside observer. The operational test is counterfactual ablation: removing the flow over the response time $\tau_d$ should break the system's self-payment loop. The postulate demarcates to which systems $\eta_L$ is applicable (bacterium, biosphere, corporate LLM stack) and to which it is not (LLM-as-agent, in whose loop self-payment is structurally absent).

Extending $N_{\text{max}}(t)$ to the integral of three flows (1) requires extending the self-payment postulate to each summand; otherwise the formula is applied ad hoc.

*Self-payment criterion for $\dot{E}_{\text{store}}, \dot{E}_{\text{grow}}$.* A flow counts as self-payment if its counterfactual removal over the response time $\tau_d$ breaks the self-payment loop — leads to a fall of $-dF_X/dt$ below threshold. This is a direct application of counterfactual ablation (Andriishin 2026, § 2.2, § 6): the $c$-test (subtracting the flow, measuring the change in the rate of free-energy depletion).

*Parameters of the criterion.* $\tau_d$ is the characteristic time of the system: bacterium $\sim 30$ min, biosphere $\sim 10$ years, LLM-corp — a quarter. The admissible telescoping of the self-payment cascade — exchange within the boundary in $\le \tau_d$.

For the biosphere and the bacterium all three flows of (1) lie within the boundary; for LLM-as-agent the self-payment loop is structurally absent (Andriishin 2026, § 3.4). The biospheric paradigm case rests on the cumulative growth of information-storage capacity through a series of evolutionary transitions with increasing carrier complexity (Maynard Smith and Szathm'ary 1995: from the RNA world to DNA–protein cells, from prokaryotes to eukaryotes, from unicellular to multicellular; in parallel (Kauffman 1993): autocatalytic networks and self-organisation as a mechanism for forming the structural $C_v^{\text{static}}$ prior to selective pressure). These transitions are the historical support for treating $\dot E_{\text{grow}}(t)$ as a self-paid component on aeonic windows.

*Caveat for the biosphere: refresh is selection, not Lemma 1.* To avoid physics-overreach, we explicitly distinguish two senses. For the bacterium $\dot M_{\text{refresh}}$ is the active holding of a bit against relaxation (Lemma 1, $k_B T\ln2$/bit as the literal payment for each cycle of rewriting a methyl mark). For the biosphere on aeonic windows $\dot M_{\text{refresh}}$ is *selection and replication* (updating of genomic content through differential reproduction), and *not* refresh in the sense of Lemma 1: the genome is not held by active correction of every bit against thermodynamic erosion over the window $\tau_{1/2}$. Lemma 1 applies to the genome only as an *order-of-magnitude bound* on the minimal cost of copying information, not as a literal refresh count. The relational structure "memory–environment" and the asymptotics of Lemma 2 carry over to the biosphere; the Landauer scale of Lemma 1 — only in the order-of-magnitude sense.

*Classical test cases for the agency-criterion from philosophy of biology.* The biological paradigm cases naturally extend to the classical borderline cases of theoretical biology, where the question of the boundary of the agent and the responsible bearer of self-payment has a long tradition of discussion. (i) *Endosymbiosis* (a mitochondrion in a eukaryotic cell): does self-payment refer to the cell as a whole or to the mitochondrion separately? The mitochondrion retains a reduced genome and its own expression apparatus, but its $\dot E_{\text{store}}$ is paid for by metabolites imported from the host — telescoping of the self-payment cascade within the cell. (ii) *Multicellularity* (cell vs organism): a single cell that has lost the capacity for autonomous reproduction within a tissue has self-payment only at the level of the organism; the counterfactual removal of the organism breaks the self-payment loop over a $\tau_d$ shorter than the counterfactual removal of the cell. (iii) *Holobiont* (host + microbiome): self-payment may be satisfied at two distinct boundaries, the $c$-test separates the cases of facultative and obligate symbioses. (iv) *Eusocial insects* (worker vs colony): a worker ant as a reproductive zero has self-payment only at the colonial boundary; the unit of selection and the unit of self-payment structurally coincide. (v) *Lichen* (fungus + photobiont): two boundaries or one? These cases are natural tests of the resolving power of the self-payment criterion in the non-stationary regime (Maynard Smith and Szathm'ary 1995; Kauffman 1993).

**LLM-corp (an extended illustrative boundary).** A non-biological borderline case: with the operating corporation included, $\dot{E}_{\text{store}}$ is paid for by revenue, and its counterfactual removal over a quarter breaks the self-payment loop. Self-payment is then satisfied at the level of the corporation, and $\eta_L^{\text{corp}}(t) > 0$ is defined — extending the boundary extends the responsible subject (Andriishin 2026, § 2.2 (the boundary-escalation principle), § 3.4 (escalating boundaries $B_1 \subset B_2 \subset B_3 \subset B_4$ and the order-of-magnitude estimate $\eta_L^{\text{corp}} \lesssim 10^{-25}$); Supplementary § S3.4 of paper #1 — full component-wise analysis).

A summary by paradigm classes (biological — central; LLM classes — illustrative borderline cases for testing the self-payment criterion on non-biological systems):

| System class | Boundary | Source of $\dot{E}$ | Self-payment | $\eta_L(t)$ |
|---|---|---|---|---|
| Bacterium in a chemostat | Cell | Internal ATP (methylation, repair, replication) | Satisfied | Defined |
| Biosphere over aeons | Biosphere | Solar exergy processed by photosynthesis within the boundary | Satisfied | Defined |
| LLM-as-agent (boundary) | Weights at inference | Data centre (external) | Structurally absent | undefined |
| LLM-corp (boundary) | LLM + operating corporation | Corporate revenue | Satisfied at the corporate level | $\lesssim 10^{-25}$ |

*Status of the numerical values in the table.* The specific estimates ($\tau_d \sim 30$ min for the bacterium, $\sim 10$ years for the biosphere as the response time of biogeochemical cycles — *not* the aeonic window of accumulation of $|M_{\le t}|$, $\sim$ a quarter for the corporation) are order-of-magnitude scales, not calibrated measurements. Full derivations of $\eta_L^{\text{corp}}$ through the capacity bound of the numerator and the corporate exergy denominator — (Andriishin 2026, § 3.4 + Supplementary § S3.4); the biospheric $\tau_d$ is the characteristic response time of the CO₂ cycle, distinct from the $\tau_E$ of climate drift and from the aeonic window of accumulation of genetic information (Landenmark et al. 2015; Maynard Smith and Szathm'ary 1995). The numbers in the table are used for a qualitative separation of the paradigm cases (defined / undefined / orders of magnitude), not for quantitative predictions.

*Relation to closure traditions.* Counterfactual ablation (the $c$-test) for self-payment is structurally close to operational closure in autopoiesis (Maturana and Varela 1980), closure of constraints (Mossio et al. 2016), and closure to efficient causation in Rosen's M,R-systems (Rosen 1991). This is a division of labour, not a displacement: the closure traditions articulate the *ontological* structure of the closedness of the self-payment loop, whereas the $c$-test adds an *operational* criterion for its diagnosis — a measurable threshold $\tau_d$ for $-dF_X/dt$, without commitments about the specific form of the causal structure (closure to efficient causation in Rosen, circular causation in autopoiesis, the ontological primacy of constraints in Mossio–Montévil–Longo). The positioning is symmetric to the deferential treatment of the same traditions in (Andriishin 2026, § 4.7), not supersessionist, and is consistent with the demarcation debate on the minimal criteria of life (Cleland 2019).

*Applicability principle for (1).* Formula (1) and $\eta_L(t)$ are applicable only under self-payment for all three flows; otherwise $\eta_L(t)$ is assigned to the extended boundary with an explicit redirection. The connection with the Pearl/Friston blanket (Bruineberg et al. 2018; Bruineberg et al. 2022) — § 7.

## 3. Predictive information under growing memory

### 3.1 Extended definition

The stationary predictive information (Andriishin 2026, formula (1a) § 2.1) is defined as

$$I_{\text{pred}}^{\text{stat}}(X, \tau) := I\!\left(M_t;\, X_E^{[t, t+\tau]}\right),$$

where $M_t$ is the configuration of the system's internal degrees of freedom at time $t$. In the stationary regime $M_t$ is a fixed object; its statistics do not change over time.

In the non-stationary regime with growing memory, $M_t$ is the current, actually updated layer, whereas the accumulated historical layers form the archival memory $M_{\lt t}$. The full informational structure of the system is formalised as a *joint random variable*

$$M_{\le t} := (M_t,\, M_{\lt t}),$$

defined on the product of the measurable spaces of the current and archival layers. The joint distribution $p(M_t, M_{\lt t})$ automatically accounts for possible correlations between the layers (in particular, when they physically overlap).

*Core definition of the non-stationary predictive information.* As the core definition of the non-stationary predictive information in the present work we take the *one-step mutual information* between the current environmental state and its state at the horizon $\tau$ conditioned on the current environment estimate $\hat P(t)$:

$$I_{\text{pred}}(t, \tau) := I\!\left(X_t;\, X_{t+\tau} \mid \hat P(t)\right). \tag{4}$$

This is an operational, step-wise characteristic of the system's predictive capacity: it is measurable in any system with a learner (biosphere, LLM, bacterium) without access to the latent environmental parameter $\theta(t)$, through a direct estimate of the mutual information between the environmental trajectories and the observer's model $\hat P(t)$. By construction $I_{\text{pred}}(t, \tau) \le \ln k$ (for a finite alphabet of $k$ states). The structural correspondence with the stationary predictive information (Andriishin 2026, (1a) § 2.1) is the one-step projection of the interval-based MI $I(M_t; X_E^{[t,t+\tau]})$; exact equality is attained for a first-order Markov chain with $M_t = X_t$ at $\tau = 1$ and $\hat P(t) \to P_0$, $\theta(t) \equiv \theta_0$ (§ 2.2, *Difference of operationalisations*).

*Connection with memory.* The estimate $\hat P(t)$ is a measurable function of $M_{\le t} = (M_t, M_{\lt t})$; by the DPI on the chain $M_{\le t} \to \hat P(t) \to (X_t, X_{t+\tau})$ we have $I(X_t; X_{t+\tau}\mid\hat P(t)) \le I(M_{\le t};\, X_E^{[t,t+\tau]})$ given the sufficiency of $\hat P(t)$ relative to $X_E^{[t,t+\tau]}$ (Cover and Thomas 2006, ch. 2.8–2.9); equality — under optimal coding of the sufficient statistic. The decomposition $M_{\le t} = (M_t, M_{\lt t})$ retains meaning: $M_t$ is the current updated layer, $M_{\lt t}$ the archive. Examples: bacterium — $M_t$ receptor methylation, $M_{\lt t}$ genotype; LLM — $M_t$ KV-cache, $M_{\lt t}$ weights; civilisation — $M_t$ administrative procedures, $M_{\lt t}$ institutional traditions.

*Categorial status of the model.* The "model of the environment" $\hat P(t)$ in the present work is a *functional kind*: it is realised in any substrate satisfying (i) a measurable mutual information $I(\hat P(t);\, X_E) > 0$ with the environmental trajectory and (ii) the self-payment criterion of § 2.6. This is an operational definition, not an ontological identification of the model with a specific substrate. A bacterium (receptor methylation), an LLM (KV-cache + weights), a civilisation (administrative procedures + institutional traditions) are *different substrates of one functional structure*, not ontologically identical systems; what they have in common is the structure "the state of an internal degree of freedom that is statistically informative about future environmental states and held against thermodynamic erosion at the expense of the boundary's own budget". This move is consistent with the functionalist tradition in the philosophy of mind/biology, where multiple realizability separates functional structure from substrate realisation (used in the sense of functional analysis (Cummins 1975), not in the sense of Putnam's machine-state functionalism; the distinction natural kind / functional kind follows (Boyd 1991; Khalidi 2013) — a functional kind is individuated by relational structure, not by microphysical homogeneity).

*Remark on the excess-loss form.* The alternative through the cumulative regret $L_{\text{excess}}(t) = \sum_{s\le t} [\ln p(X_s\mid X_{s-1},\theta(s)) - \ln \hat p(X_s\mid X_{s-1},\hat\theta(s))]$ is a kindred BNT Class II object (Bialek et al. 2001), operational only in simulations with ground truth $\theta(s)$; in real systems $\theta(s)$ is inaccessible. The connection with (4), the asymptotics, and the adiabatic limit — Supplementary § S8.

### 3.2 Separation into current and archival layers

The chain rule of MI gives an exact decomposition of the information about the environment's future contained in the full memory:

$$I\!\left(M_{\le t};\, X_E^{[t, t+\tau]}\right) = I\!\left(M_t;\, X_E^{[t, t+\tau]}\right) + I\!\left(M_{\lt t};\, X_E^{[t, t+\tau]} \mid M_t\right). \tag{5}$$

The first summand is the "fresh" predictive information of the current layer (for $M_{\lt t} = \emptyset$ — the original predictive information Andriishin 2026); the second is the conditional contribution of the archival layer, nonnegative and vanishing when the archive is redundant. By the DPI: $I_{\text{pred}}(t,\tau) \le I(M_{\le t};\, X_E^{[t, t+\tau]})$, with equality attained when $\hat P(t)$ is sufficient as a function of $M_{\le t}$. The decomposition (5) is the basis for the formulation of Lemma 2 § 4.4, which separates the contributions of the updated and the ageing fractions of memory.

*Upper bound on predictive information.* For a finite state alphabet $|\Omega| = k$ and any estimate $\hat P(t)$, we have $I_{\text{pred}}(t, \tau) \le \ln k$ — the operational one-step MI is bounded above by the entropy of the chain's states. Separately: for an ergodic chain the stationary limit of the one-step $I_{\text{pred}}^{\text{stat}}(\tau)$ is majorised by the excess entropy $C_\mu$ (Crutchfield and Feldman 2003), which may exceed $\ln k$ for processes with memory greater than first order; the proper bound $\ln k$ applies precisely to the one-step conditional MI.

### 3.3 Reducibility to the stationary definition

In the stationary Markov regime $\theta(t) \equiv \theta_0$, the estimate $\hat P(t) \to P_0$ over a time of order $\tau_{\text{relax}}$, and the one-step MI (4) reaches a finite constant $I_{\text{pred}}^{\text{stat}}(\tau) = I(X_t; X_{t+\tau} \mid \theta_0)$ — the stationary value, consistent with paper #1 in the sense clarified in § 2.2 (the one-step projection of the interval-based $I(M_t; X_E^{[t,t+\tau]})$ Andriishin 2026, (1a) § 2.1). The non-triviality of (4) manifests in the non-stationary regime: in real systems optimal forgetting is not achieved (bacterial methylation lags behind the ligand; LLM weights are not updated for vocabulary drift; the archive of institutional norms is redundant), and $M_{\lt t}$ has a nonzero but non-predictive contribution to (5). This prepares the definition of nostalgia (§ 4): the archival layer remains informationally rich but grows poorer in predictive value owing to the environmental shift.

## 4. The phenomenon of nostalgia

### 4.1 Definition

*Ontological status.* Nostalgia is a *relational property of the "memory–environment" pair over time*: a mismatch between the statistics held by the system and the optimal predictive statistics relative to the current environmental state. It is not a property of memory as an isolated object (the same fixed set of bits in the same substrate may be fully predictive relative to one environment and fully nostalgic relative to another) and not a property of the environment (the same environmental drift produces different nostalgia in systems with different refresh rates). Analogously to the distinction between the *referent* and the *operationalisation* of the vital scale (Andriishin 2026, § 1, § 2.7), the referent of nostalgia (a relational property) differs from its two operationalisations: $\nu^{\text{theor}}$ through the counterfactual oracle $\theta(t)$ (a theoretical instrument) and $\nu^{\text{op}}$ through the empirical excess entropy (a measurable procedure). Both operationalisations are estimates of one referent; the DPI inequality $\nu^{\text{op}} \le \nu^{\text{theor}}$ is a structural relation between two estimates, not a definition of the object.

*Categorial status and multiple realizability of nostalgia* (parallel to § 3.1 *Categorial status of the model*). Nostalgia as a relational property of the "memory–environment" pair is a functional kind, multiply realisable: it is realised in any pair satisfying (i) a measurable memory refresh rate $\dot{M}_{\text{refresh}} \lt \infty$; (ii) a measurable environmental drift with finite $\tau_E$; (iii) the self-payment criterion of § 2.6 for holding memory. Biological (bacterial methylation — ligand gradient), computational (LLM-corp weights — query statistics), and social (institutional norms — socio-economic reality) realisations are *different substrates of one relational structure*. This does **not** assert the ontological identification of the informational nostalgia of the present work with psychological nostalgia (the everyday word is retained for intuitive support, but the referent is strictly distinct: psychological nostalgia is a phenomenal state of a subject, informational nostalgia is a mismatch of statistics; see § 8.2 item (ii)).

*Terminological risk.* The choice of the word "nostalgia" over neutral alternatives (memory staleness, predictive obsolescence) is a methodological risk of smuggling everyday intuition into the formalism, analogous to the risk diagnosed (Bruineberg et al. 2022) for the term "Markov blanket". The word is retained for intuitive support (the biological recognisability of the phenomenon) and continuity with the population-biology tradition of bet-hedging (Kussell and Leibler 2005). The formalism (7)/(7') does not depend on the choice of word; a reader who sees a risk in this choice may replace the term with predictive obsolescence without loss of content. The everyday psychological meaning is a categorially distinct referent.

Nostalgia $\nu(t)$ is an *operational measure* of the said relational property: the fraction of held memory that no longer reduces the discrepancy between the model and the environment (*the fraction of potentially available information about the current environment that is not used for prediction*). A *theoretical* and an *operational* normalisation are distinguished through the one-step (4); the separation is essential for consistency with the protocol of § 5.2 / Supplementary § S5.2 item 5.

*Theoretical nostalgia* $\nu^{\text{theor}}(t)$ — through an ideal observer with full sufficient statistics relative to the true latent $\theta(t)$:

$$\nu^{\text{theor}}(t) := 1 - \frac{I_{\text{pred}}(t, \tau)}{I_{\text{pred}}^{\text{opt}}(t, \tau)}, \quad I_{\text{pred}}^{\text{opt}}(t,\tau) := I\!\left(X_t;\, X_{t+\tau} \mid \theta(t)\right). \tag{7}$$

The denominator $I_{\text{pred}}^{\text{opt}}$ requires access to $\theta(t)$ and is therefore *counterfactual* (like $L_{\text{excess}}$ in Supplementary § S8): computable in simulations with ground truth, not directly measurable in real systems. $\nu^{\text{theor}}$ is the measure for rigorous statements (Lemma 2, $\nu_c$, $C_v^{\text{predictive}}$).

*Operational nostalgia* $\nu^{\text{op}}(t)$ — through an empirical upper bound without the oracle $\theta(t)$:

$$\nu^{\text{op}}(t) := 1 - \frac{I_{\text{pred}}(t, \tau)}{C_\mu^{\text{emp}}(t, \tau)}, \tag{7'}$$

where $C_\mu^{\text{emp}}(t, \tau)$ is the empirical excess entropy from the trajectory over a sliding window $\tau_w \gg \tau_{\text{relax}}$ (Crutchfield and Feldman 2003); for a slowly drifting environment with $\tau_w \ll \tau_E$ it gives an unbiased estimate of $I_{\text{pred}}^{\text{opt}}$. $\nu^{\text{op}}$ is the figure in the operational predictions (P.1, P.2 in § 8.6), accessible to measurement in any system with a learner.

For both normalisations $\nu \in [0, 1]$ — the DPI on the chain $\theta(t) \to \hat P(t) \to (X_t, X_{t+\tau})$: conditioning on the coarse-grained $\hat P(t)$ does not give MI greater than conditioning on the sufficient one (Cover and Thomas 2006, ch. 2.8–2.9). The poles: $\nu = 0$ — the estimate reaches the power of the optimal observer; $\nu = 1$ — no contribution beyond the unconditional one. The denominator — the one-step conditional MI — is bounded above by $\ln k$ (Cover and Thomas 2006).

*Convention on the boundary of the $\nu$ domain.* By analogy with the convention on the boundary of the $\eta_L$ domain in paper #1 (Andriishin 2026, § 2.1), three boundary regimes are fixed in which (7)/(7') do not define a scale: (i) for $I_{\text{pred}}^{\text{opt}}(t, \tau) = 0$ we set $\nu(t) := 0$ as the limit $\lim_{I^{\text{opt}} \to 0^+} \nu$ with simultaneous $I_{\text{pred}} \to 0$; (ii) for $I_{\text{pred}}^{\text{opt}}(t, \tau) \le \varepsilon_{\text{floor}}$, where $\varepsilon_{\text{floor}}$ is the *floor* of the MI estimator (for the KSG estimate a typical threshold is $\sim 0.01$ nat at limited sample length Kraskov et al. 2004), $\nu(t)$ is set *undefined*, and Lemma 2 is inapplicable in this regime (the noise/signal ratio in (7) exceeds the controlled threshold); (iii) for $|M_{\le t}| < N_{\min}$ — on short initial intervals $t < \tau_d$, where the size of the active memory has not reached the power-analysis threshold of the protocol (Supplementary § S5.2) item 1 — the operational nostalgia $\nu^{\text{op}}(t)$ is set undefined, and the main-text $\nu(t)$ refers to the asymptotic regime $t \gg \tau_d$ with $|M_{\le t}| \ge N_{\min}$. This discipline is symmetric to the convention of paper #1 on removing the system from the domain of the scale in the absence of paid dissipation.

*When the index is omitted.* $\nu(t)$ is $\nu^{\text{theor}}$ in theoretical statements (Lemma 2, $\nu_c$, $C_v^{\text{predictive}}$) and $\nu^{\text{op}}$ in operational predictions and the illustrations of § 6 (where the ground truth of the simulation makes both normalisations coincide).

*Stability of $\nu_c$ under a change of normalisation.* By the DPI (see § 3.2) $C_\mu^{\text{emp}}(t,\tau) \le I_{\text{pred}}^{\text{opt}}(t,\tau)$ given a sufficient estimation window; a smaller denominator in (7') entails a larger fraction $I_{\text{pred}}/C_\mu^{\text{emp}}$, hence $\nu^{\text{op}}(t) \le \nu^{\text{theor}}(t)$ pointwise. In this work a *working assumption* is adopted: $\nu_c^{\text{theor}}$ as the regime threshold (§ 5.3) is stable under the transition to $\nu^{\text{op}}$, with a possible numerical shift $\nu_c^{\text{op}} \ne \nu_c^{\text{theor}}$. A rigorous stability theorem (monotonicity of the regime threshold under the DPI inequality) is an open problem § 8.5; the pre-registered protocol § 8.6 (P. 1) is formulated in terms of $\nu_c^{\text{op}}$ as an operationally measured threshold.

*Terminology.* "Nostalgia" is a technical term denoting a relational property of the "memory–environment" pair (see *Ontological status* above); the measures (7)/(7') are its operationalisations. In EN — *informational nostalgia*.

*Connection with virtual concept drift.* $\nu(t)$ is a thermodynamic operationalisation of virtual concept drift from ML (Schlimmer and Granger 1986; Tsymbal 2004; Gama et al. 2014): $X_E$ drifts, the optimal $M^{\text{opt}}$ shifts, while the held $M_{\le t}$ remains fixed. ML concept-drift detectors (ADWIN (Bifet and Gavald`a 2007), DDM Gama et al. 2004) rely on statistical tests; the present work adds a thermodynamic scale through (7)–(8) and relates the detection rate to the $\dot{M}_{\text{refresh}}$ of Lemma 2.

### 4.2 Thermodynamic cost

Nostalgic bits are bits paid for by $E_{\text{store}}$ (by Lemma 1) and contributing nothing to $I_{\text{pred}}$. If $|M_{\le t}|$ is the full memory size, the nostalgic part has size $\nu(t) \cdot |M_{\le t}|$ and holding it against degradation requires power

$$\dot{E}_{\text{store}}^{\text{nostalg}}(t) \ge \frac{\nu^{\text{op}}(t) \cdot |M_{\le t}|}{\tau_{1/2}} \cdot k_B T \ln 2. \tag{8}$$

(The context is operational — the measurable payment for holding the non-predictive fraction of memory, so (8) uses $\nu^{\text{op}}$ by (7').) *Consistency with (2).* Formula (8) is a coarse-grained bit-level estimate in which the system accuracy threshold $\varepsilon$ is fixed by the architecture of the memory carrier (typical biological memory operates at $\varepsilon \sim 1/e$, giving a multiplicative factor of order unity). Upon explicitly unfolding the dependence on $\varepsilon$ through (2), the right-hand side of (8) is multiplied by $\ln 2/(2\varepsilon)$; for systems with variable $\varepsilon$ this refinement is essential. The summand enters the total $E_{\text{store}}$ (1) — the denominator of $\eta_L(t)$ — without contributing to the numerator: nostalgia monotonically lowers $\eta_L$. *Invariance in $\varepsilon$.* The qualitative conclusion $\eta_L \to 0$ does not depend on the choice of refresh scheme — naive refresh (2) or majority-vote (2'): both bounds are linear in $t/\tau_{1/2}$, so $\varepsilon$ affects only the prefactor of holding, not the efficiency collapse itself.

*Remark (Sagawa–Ueda and the self-payment loop).* A nostalgic bit is correlated with the past environment ($C_\mu > 0$), so formally the holding of correlated memory could be made cheaper, à la Sagawa–Ueda, relative to the naive-refresh bound. In the self-payment loop there is no such loophole: the acquisition and holding of this correlation have already been paid for by the same internal budget (by the measurement–erasure symmetry, the payment $\ge k_B T \cdot I$ for the correlation itself), and crediting a Sagawa–Ueda mitigation on top of this is double counting, equivalent to a Maxwell demon inside the system's boundary. Therefore (8) is a naive-refresh floor (as is Lemma 1), not an estimate with feedback-cheapening already accounted for.

*Connection with Still et al.* (8) is a non-stationary analogue of the Still–Sivak–Bell–Crooks 2012 inequality in instantaneous form $\dot W_{\text{diss}} \ge k_B T \ln 2 \cdot (\dot I_{\text{mem}} - \dot I_{\text{pred}})$ (bit notation, see Still et al. 2012, eq. 3–6); in our case it operates in the limit $\dot I_{\text{pred}}^{\text{contrib}} \to 0$ of the nostalgic layer — ageing bits contribute to $\dot I_{\text{mem}}$ (they are maintained) but not to $\dot I_{\text{pred}}$, which is fixed by the unfolding through $\tau_{1/2}$ in (8). Information reservoirs — (Mandal and Jarzynski 2012; Boyd et al. 2017); bipartite TE — (Hartich et al. 2014).

### 4.3 Methodological refinement of $C_v$

Distinguishing the currently predictive and the nostalgic parts of memory requires a refinement of the operationalisation of $C_v$ (Andriishin 2026, § 2.3). In the stationary regime two classes of proxies — entropy rates $h_\mu$ (Lempel–Ziv, permutation entropy) and predictive information (excess entropy, $C_\mu$ Crutchfield and Feldman 2003) — give consistent characteristics. In the non-stationary regime the distinction becomes diagnostic:

- $C_v^{\text{static}}$ — *structural complexity* through proxies of the first class (Lempel–Ziv, assembly index
  Sharma et al. 2023); does not distinguish the predictive and the nostalgic parts.
- $C_v^{\text{predictive}}$ — *dynamic predictive power* through proxies of the second class, computed
  relative to the current environmental statistics.

Connection:

$$C_v^{\text{predictive}}(t) = (1 - \nu(t)) \cdot C_v^{\text{static}}(t) \tag{9}$$

to first approximation (uniform distribution of predictiveness across bits). The refinement (9) does not annul the triad $T_v, C_v, S_v$: in the stationary limit $\nu \to 0$ the distinction disappears. In the non-stationary regime, the divergence of $C_v^{\text{static}}$ and $C_v^{\text{predictive}}$ is a diagnostic signal of nostalgia.

### 4.4 Lemma on the inevitability of nostalgia

The central result of the section is a conditional lower bound for the asymptotic fraction of nostalgia in the slow-driving class, with an explicit mechanism of emergence through the desynchronisation of the optimal predictive statistics. In the general non-stationary case a universal lower bound is underivable; the lemma is formulated conditionally for the class of § 5.2.

*Precedent.* Lemma 2 is a non-stationary analogue of Still–Sivak–Bell–Crooks (Still et al. 2012) with an explicit unfolding through the memory refresh rate; a parallel is the stochastic thermodynamics of learning (Goldt and Seifert 2017). The extension relative to (Still et al. 2012): (i) the tie to the characteristic drift time $\tau_E$ and $\dot{M}_{\text{refresh}}$ through $(*)$; (ii) the decomposition of $M_{\le t}$ into current and archival layers (§ 3.1); (iii) the asymptotics $\liminf \nu(t)$ through $c, \tau_E$.

**Lemma 2 (inevitability of nostalgia in the OU class of slow drift).** *Let the environment $X_E$ belong to the class of § 5.2 (slow continuous drift of the logits with $\tau_E = 1/\lambda$, slow-driving $\lambda \tau_{\text{relax}} \ll 1$); let the memory be updated at a rate $\dot{M}_{\text{refresh}} \lt \infty$ by FIFO (the age distribution is asymptotically uniform for $t \ge t_0 + |M_{\le t}|/\dot{M}_{\text{refresh}}$); $I_{\text{pred}}(t,\tau)$ given by (4), $\nu^{\text{theor}}(t)$ by the normalisation (7); the per-bit-uniformity assumption of Remark 5 holds (additivity of MI over independently updated bits with constant $1+o(1)$). Introduce the computable constant*

$$c := \frac{\dot{M}_{\text{refresh}} \cdot \tau_E}{|M_{\le t}|} \tag{*}$$

*and assume $c \lt 1$ for all $t > T_0$ (otherwise Lemma 2 is inapplicable: the "sliding window" regime, see Remark 2). Then asymptotically*

$$\liminf_{t \to \infty} \nu^{\text{theor}}(t) \ge 1 - c. \tag{8a}$$

*What is non-trivial in Lemma 2.* Formally the arithmetic $c = \dot{M}_{\text{refresh}} \cdot \tau_E / |M_{\le t}|$ from $(*)$ is simply the fraction of refresh updates over the environment's correlation interval; the statement $\nu \ge 1 - c$ **does not follow** from this definition without two additional steps.

(i) **The OU spectral gap** (Pavliotis 2014, ch. 4) (an alternative path through MI tensorisation over the independent OU coordinates, Supplementary § S1.5) ensures the vanishing of the conditional MI of the stale fraction $I(M_{\le t}^{\text{stale}};\, X_E^{[t,t+\tau]} \mid M_{\le t}^{\text{refr}}) \to 0$ at an exponential rate $e^{-2\Delta t/\tau_E}$ for $\Delta t \gg \tau_E$ (the logit covariance decays as $1/\tau_E$, but the conditional MI as $2/\tau_E$, since $I = \tfrac12\rho^2$ for a Gaussian pair, $\rho = e^{-\Delta t/\tau_E}$; rigorous derivation — Supplementary § S1.5). This is a substantive *probability-theoretic* step: the decorrelation of the OU semigroup between the moment of fixing the stale fraction and the current environmental state is a non-trivial property of the OU class, *not* valid in general for slow-driving processes with long memory (1/f noise, fractional OU, multi-scale drift).

(ii) **Per-bit-uniformity** (Remark 5, Supplementary § S5.1) ensures the additivity of the mutual information over coordinates: $I_{\text{joint}} = \sum_k I_{\text{coord}}^{(k)}$ for $K$ independently updated OU logits. Without per-bit-uniformity the total MI is bounded above by $K \cdot I_{\text{coord}}^{\max}$, which may be considerably smaller than the sum; for the transition from the arithmetic fraction $c$ to the information-theoretic inequality $\nu \ge 1-c$ exactly this additivity is required.

Without both steps (i)+(ii), Lemma 2 is a definition disguised as a theorem. Since both premises hold *precisely for the OU class* (Remarks 4–6 § S5.1) and are *not proved in general for the slow-drift class* (Remark 4 explicitly notes: the PSP surrogate and multi-scale drift reproduce the result illustratively, not rigorously), the formal domain of applicability of Lemma 2 is the OU class of § 5.2; the carry-over to other slow-driving classes is an open problem § 8.5.

**Proof.** The principal condition is the *exponential decorrelation of the optimal predictive statistics with characteristic time $\tau_E$* in the slow-drift class. The optimal predictive feature $\xi(t)$ is a projection of $M_{\le t}$ that is statistically sufficient relative to $X_E^{[t, t+\tau]}$. In the OU class of § 5.2 the autocorrelation $\mathbb{E}[\theta(t)\theta(t-s)] \propto e^{-s/\tau_E}$ is a standard result for OU (Remark 6 fixes the transfer of the decay through softmax to $\xi$). By the DPI, a feature without updating over a horizon $\gg \tau_E$ belongs to the nostalgic fraction by (7).

Over an interval $\tau_E$ the system updates no more than $\dot{M}_{\text{refresh}} \cdot \tau_E$ bits; under FIFO, for $t \ge t_0 + |M_{\le t}|/\dot{M}_{\text{refresh}}$ the age distribution is asymptotically uniform, and the fraction of un-updated bits over the interval $\tau_E$ is bounded below by $1 - c$ from $(*)$. For an adversarial schedule — $1 - c'$ with $c' \in (c, 1)$; (8a) is formulated for FIFO.

In MI terms, the decomposition (5) is applied to the partition $M_{\le t} = (M_{\le t}^{\text{refr}}, M_{\le t}^{\text{stale}})$ by the moments of last update:

$$I(M_{\le t};\, X_E^{[t,t+\tau]}) = I(M_{\le t}^{\text{refr}};\, X_E^{[t,t+\tau]}) + I(M_{\le t}^{\text{stale}};\, X_E^{[t,t+\tau]} \mid M_{\le t}^{\text{refr}}),$$

where the first term is the MI of the updated fraction, the second the conditional MI of the stale fraction (chain rule). The key step is $I(M_{\le t}^{\text{stale}};\, X_E^{[t,t+\tau]} \mid M_{\le t}^{\text{refr}}) \to 0$ for $\Delta t \gg \tau_E$. $M_{\le t}^{\text{stale}}$ was fixed at a moment $s \le t - \Delta t$ on the basis of $X_E^{[\,\cdot\,, s]}$ at $\theta(s)$; by the Markovianity of the pair $(\theta, X_E)$ the chain $M_{\le t}^{\text{stale}} \to (\theta(s), X_E^{[\,\cdot\,,s]}) \to \theta(t) \to X_E^{[t,t+\tau]}$ is a correct DPI structure, whence
$$I\bigl(M_{\le t}^{\text{stale}};\, X_E^{[t,t+\tau]} \,\bigm|\, M_{\le t}^{\text{refr}}\bigr) \le
 I\bigl((\theta(s), X_E^{[\,\cdot\,,s]});\, \theta(t) \,\bigm|\, M_{\le t}^{\text{refr}}\bigr).$$
In OU $\mathbb{E}[\theta(t)\theta(s)] \propto e^{-(t-s)/\tau_E}$; the dependence of $X_E^{[\,\cdot\,,s]}$ on $\theta(t)$ is transmitted only through $\theta(s)$. By the standard result on the spectral gap of the OU semigroup (the Bakry–Émery functional inequality for a Gaussian measure; see Pavliotis 2014, ch. 4, the theorem on the spectral decomposition of OU and exponential ergodicity) the right-hand side $\le C\,e^{-2\Delta t/\tau_E}$ (the MI decays as $2/\tau_E$, twice as fast as the covariance of $\theta$; Supplementary § S1.5). Averaging over ages under the uniform FIFO distribution gives the vanishing of the contribution of the stale fraction in (5) to accuracy $O(c)$ for $c \ll 1$; substitution into (7) under the assumption of Remark 5 gives (8a). The chain rule on MI is a standard identity; transfer entropy (Schreiber 2000) is not invoked. $\square$

**Remarks (brief overview).** Six Remarks to Lemma 2 unfold the conditions of its applicability: (1) the stationary limit is set by $\dot{E}_{\text{store}} = \dot{E}_{\text{grow}} = 0$, not by $\tau_E \to \infty$; (2) the constant $c$ from $(*)$ is computable through the measurable $(\dot{M}_{\text{refresh}}, \tau_E, |M_{\le t}|)$, at $c \ge 1$ — the sliding-window regime; (3) (8a) is substantive beyond the arithmetic definition of $c$ under the conditions of decorrelation, an ideal observer, and FIFO refresh; (4) the basic Lemma 2 is proved for the OU class, the generalisation to the general non-stationary class is an open problem § 8.5; (5) per-bit-uniformity ensures the transition from the arithmetic fraction $c$ to the information-theoretic inequality $\nu \ge 1-c$; (6) softmax-regularity for OU ensures the transfer of OU decorrelation onto any bounded $\xi(\theta)$ with the same $\tau_E$. Full formulations of Remarks 1–6 with detailed arguments and counterarguments — Supplementary § S5.1.

Lemma 2 is constructive: under $(*)$ it indicates an explicit mechanism of nostalgia (desynchronisation of the optimal feature) and a physical obstacle to its being zero.

*Biological predecessor in population biology.* Lemma 2 has a direct structural predecessor in the work (Kussell and Leibler 2005), which proved, in a population-biology setting, an information inequality for bet-hedging strategies in fluctuating environments: the optimal population growth rate is related to the mutual information between the environmental signal and the phenotypic distribution, and a memory of the environment exceeding its correlation scale monotonically lowers fitness. The extension of the present work relative to (Kussell and Leibler 2005): (i) the thermodynamic scale — the transition from relative fitness to the Landauer cost of holding (8); (ii) the explicit decomposition of $M_{\le t}$ into the updated and archival layers (5); (iii) the asymptotics $\liminf \nu \ge 1 - c$ with a computable constant $c$ through $(*)$, not derivable from the bet-hedging formalism (Kussell and Leibler 2005). The structural correspondence confirms that Lemma 2 is a non-stationary thermodynamic reflex of a known population-biology regularity, not an isolated construction.

*Level mapping bet-hedging → Lemma 2.* In the bet-hedging formalism of Kussell and Leibler (2005) information is held at the population level through the phenotype distribution; in Lemma 2 — at the individual level through the memory $M_{\le t}$. The structural parallelism: the optimal switching rate of Kussell and Leibler (2005) correlates with $\dot{M}_{\text{refresh}}$, the environmental correlation scale of Kussell and Leibler (2005) with $\tau_E$, the fitness cost of stale information with $\dot{E}_{\text{store}}^{\text{nostalg}}$. This map is an illustrative analogy of two levels, not a reduction of one to the other.

*Carry-over to non-biological paradigm cases.* The carry-over of Lemma 2 to LLM-corp, cultural, and institutional memory is correct at the level of the functional structure (mismatch of the optimal statistics) but *not* at the level of the thermodynamic scale of Lemma 1 — for LLM-corp the scale is recomputed through § 2.6 (corporate exergy as a self-paid flow), for cultural memory $\tau_{1/2}$ is operationally undefined (the carrier is ambiguous: book, digital archive, neural network; see § 8.5). To non-biological paradigm cases there carry over the relational structure "memory–environment" and the asymptotics $\nu \ge 1-c$ under satisfied $(*)$; the Landauer cost of holding (8) retains meaning only under a satisfied self-payment criterion.

*Periodic reset of environment-relevant information.* Biological systems implement periodic reset of the nostalgic layer for the sake of restoring $\nu \to 0$ through several distinct mechanisms: (i) sporulation in *Bacillus subtilis* — an extreme periodic memory wipe upon environmental deterioration through the formation of an endospore with a sharp reduction of metabolic activity and a phenotype switch (Errington 2003); (ii) diapause egg formation in *Daphnia* — seasonal phenotypic plasticity with conservation of genetic material through resting stages restorable upon the return of favourable conditions (Decaestecker et al. 2007); (iii) epigenetic reprogramming in mammals — almost complete erasure of DNA methylation with re-establishment in early embryogenesis and the germ line (Reik 2007). In all three cases, in the light of Lemma 2 — periodic reset of the archival layer as a mechanism for restoring the adaptive regime in a new generation or life cycle.

## 5. The non-stationary $\eta_L(t)$ and its bounds

### 5.1 Definition and three regimes

The non-stationary Landauer efficiency is defined as

$$\eta_L(t) = \frac{I_{\text{pred}}(t, \tau)}{N_{\text{max}}(t)}, \tag{10}$$

where the numerator is given by formula (4) — the one-step predictive information under the current environment estimate $\hat P(t)$ — and the denominator by formula (1). This is an operational quantity: measurable in any system with a learner without access to the latent environmental parameter. Under (a) irreversibility, (b) accounting for all three summands of (1), (c) correct operationalisation of $I_{\text{pred}}(t,\tau)$ by (4), (d) $|M_{\le t}| \ge N_{\min}$ from the protocol item 1 (see Supplementary § S5.2), we have $\eta_L(t) \in [0, 1]$ in the asymptotic regime $t \gg \tau_d$. On short initial intervals $t < \tau_d$ with $|M_{\le t}| < N_{\min}$, $\eta_L(t)$ is set *undefined* (the convention on the boundary of the domain, symmetric to paper #1 Andriishin 2026, § 2.1).

An alternative to (10) — the instantaneous efficiency through a sliding window $\tau_w$,

$$\eta_L^{(\tau_w)}(t) = \frac{I_{\text{pred}}(t, \tau)}{N_{\text{max}}^{(\tau_w)}(t)}, \quad
N_{\text{max}}^{(\tau_w)}(t) = \frac{1}{k_B T \ln 2} \int_{t - \tau_w}^{t}\!
\left[\dot{E}_{\text{actual}}^{\text{curr}}(t') + \dot{E}_{\text{store}}(t') + \dot{E}_{\text{grow}}(t')\right] dt'. \tag{10a}$$

A local efficiency, convenient for diagnosing the steady-state regime: (10) is monotonically non-decreasing in the denominator with $t$, (10a) "forgets" the distant past and may grow after a reset.

In the non-stationary regime $\eta_L(t)$ exhibits three qualitatively distinct asymptotics depending on the relation of $\dot M_{\text{refresh}}$, $\tau_E$, and $\tau_{\text{reset}}$:

- *Stationary regime:* $\dot{M} = 0$, $\nu \to 0$; $\eta_L(t) \to \eta_L^{\text{stat}} \in (0, 1]$, structurally
  corresponding to $\eta_L^{\text{stat}}$ (Andriishin 2026, § 2.1) in the sense of § 2.2 (exact coincidence in the case
  of a first-order Markov chain with $M_t = X_t$, $\hat P \to P_0$, $\tau = \tau_d$; in the general case —
  structural correspondence through different operationalisations of the numerator). For the bacterium — the exergetic
  $\eta_L \sim 3 \cdot 10^{-8}$ or $\eta_L^{\text{comp}} \sim 3 \cdot 10^{-5}$ (§ 3.5 of the cited work).
- *Drift without reset (nostalgic collapse):* by Lemma 2, $\nu(t) > 0$ and grows monotonically; by (8)
  the nostalgic part builds up $\dot{E}_{\text{store}}$ without contributing to $I_{\text{pred}}$; $\eta_L(t) \to 0$.
  The regime of an overfitted LLM without continual pre-training, of a hypothetical biosphere without epigenetic reset.
- *Drift with reset:* memory is reset at interval $\tau_{\text{reset}}$; within a cycle $\eta_L(t)$ grows after
  reset, reaches a maximum, and decreases as nostalgia grows. $\langle \eta_L \rangle$ is stably positive
  and depends on $\tau_{\text{reset}}$ and $\tau_E$. Analogies: generational turnover with epigenetic
  reprogramming (Reik 2007); periodic fine-tuning of ML systems.

The stationary regime is a special case of (1) at $\dot{E}_{\text{store}} = \dot{E}_{\text{grow}} = 0$ (§ 2.2); in this limit Lemma 2 is inapplicable, and $\eta_L(t)$ reduces to a constant.

### 5.2 The slow-drift class: OU parameterisation and applicability conditions

As the canonical class for the numerical illustrations of § 6 and the formulation of the conjecture of Supplementary § S8.1 we adopt the *OU parameterisation of the logits:* the matrix $P(t)$ is set through the logits $\theta_{ij}(t)$ of dimension $K = k(k-1)$, satisfying the OU equation $d\theta_{ij} = -\lambda(\theta_{ij} - \theta_{ij}^*)\,dt + \sigma\,dW_{ij}$ with characteristic drift time $\tau_E = 1/\lambda$ and amplitude $\sigma$. The reconstruction of $P(t)$ is a row-wise softmax; on each slice $X_E$ is a stationary Markov process with relaxation time $\tau_{\text{relax}}$. The OU class realises a slowly varying latent parameter in the Class II setting of Bialek–Nemenman–Tishby (Bialek et al. 2001, § VI); the PSP surrogate of § 6.1 — a piecewise-stationary process with Poisson switches — is used as a simulation surrogate convenient for illustrating the "drift with reset" regime.

The general statements of the work (Lemmas 1–2, the three regimes of § 5.1, the existence of $\nu_c$ below, the methodological refinement of $C_v$ § 4.3) are formulated under the following applicability conditions: (i) the irreversible computational regime (Bennett 1973; Bennett 1982); (ii) finite $\tau_{1/2}$ of the memory carrier; (iii) finite refresh rate $\dot{M}_{\text{refresh}} < \infty$; (iv) slow-driving $\lambda \tau_{\text{relax}} \ll 1$. The OU-concentration condition ($\sigma^2/\lambda \ll 1$) is additional and is needed only for the conjecture of Supplementary § S8.1.

A detailed unfolding of the adiabatic regime under two distinct conditions ($\lambda \tau_{\text{relax}} \ll 1$ for Lemma 2 vs. $\sigma^2/\lambda \ll 1$ for the BNT Class II asymptotics of Supplementary § S8), the full formulation of Remarks 1–6 to Lemma 2, and the methodological protocol against ad hoc rescues with the justification of item 5 (the ban on revising the definition of the object of measurement) — in Supplementary § S5.

### 5.3 Existence of $\nu_c$ and the optimal $\tau_{\text{reset}}^*$

$\nu_c^{\text{theor}}$ is the critical level of nostalgia (on the theoretical scale (7)) at which the system switches from the adaptive regime into the collapse regime. The balance estimate from the condition of equality of the gain $\dot{I}_{\text{pred}}$ and the losses:

$$\nu_c^{\text{theor}} \cdot |M_{\le t}| / \tau_{1/2} = (1 - \nu_c^{\text{theor}}) \cdot \dot{M}_{\text{refresh}}. \tag{11b}$$

(11b) uses two scales: $\tau_E$ (in $(*)$) and $\tau_{1/2}$ (in Lemma 1). Direct substitution in the simplest form $\nu_c^{\text{theor}} = \eta_M/(1+\eta_M)$ with $\eta_M := \dot{M}_{\text{refresh}} \tau_{1/2} / |M_{\le t}|$ gives values inconsistent with § 6.1; (11b) is treated as a *structural balance condition*, not a quantitative predictor — the work postulates the existence of $\nu_c^{\text{theor}} \in (0, 1)$ as a fixed point. The empirical proxy $\nu_c^{\text{emp}}$ of § 6.1 is the $\nu_c^{\text{op}}$ by (7'); the derivation of a quantitative formula $\nu_c^{\text{theor}}(\lambda, \sigma, K, \tau_{1/2}, \dot M_{\text{refresh}})$ — § 8.5.

The optimal $\tau_{\text{reset}}^*$ maximises $\langle \eta_L \rangle$: frequent reset increases the payment, rare reset leads to collapse:

$$\tau_{\text{reset}}^* \sim \tau_E \cdot \log(|M_{\le t}^{\text{cap}}|/\dot{M}_{\text{refresh}} \tau_E), \tag{12}$$

where $|M_{\le t}^{\text{cap}}|$ is the structural capacity limit. The biological analogy is the generation time of a taxon (lifespan correlates with the rate of change of the ecological niche).

A detailed justification of the structural status of (11b), the analogy with spin glasses (Hopfield 1982; Amit et al. 1985; Engel and Van den Broeck 2001), the connection of $\tau_{\text{reset}}^*$ with dynamic regret (Besbes et al. 2015) and no-free-lunch (Wolpert 1996) — Supplementary § S5.

## 6. Numerical illustration of the theory on Markov chains with drift

*Convention of § 6.* In the simulations of this section the ground truth $\theta(t)$ is known by construction, so $\nu^{\text{theor}}(t) = \nu^{\text{op}}(t)$; $\nu(t)$ without an index is used. The adiabatic auxiliary series (a parametric OU scan in the OU-concentration regime $\sigma^2/\lambda \ll 1$ with direct measurement of the cumulative excess-loss $L_{\text{excess}}(t)$ relative to the oracle) belongs to the open conjecture of Supplementary § S8.1 and is deferred to Supplementary § S8.3.

### 6.1 PSP simulation: Markov chains with Poisson drift

The PSP simulation (a piecewise-stationary process, $k = 8$, $\alpha = 0.3$, $\lambda = 10^{-3}$, $\tau_E = 10^3$ steps) is a numerical illustration of the three regimes of behaviour of $\eta_L(t)$ (stationary/collapse/reset) from § 5.1, as well as of the structural statement of Lemma 2 about $\nu(t) \to 1$ in the collapse regime. The class by construction satisfies the slow-driving condition ($\lambda \tau_{\text{relax}} \approx 3 \cdot 10^{-3} \ll 1$, see § 5.2); the OU-concentration condition for PSP is undefined (Poisson switches have no continuous $\sigma$), so the conjecture of Supplementary § S8.1 is not tested in PSP. The observer's memory is $\hat P(t)$ over a sliding window $\tau_w$ with Laplace smoothing, $|M_t| = k^2 = 64$.

*Fig. 1: $\eta_L(t)$ for the three regimes* (`paper/figs/Fig1.pdf`). The stationary regime reaches a horizontal asymptote $\eta_L \approx 6.1 \cdot 10^{-5}$ ($\lambda = 0$, $\tau_w = 100$, $T = 10^4$). These values are small on the absolute scale but not anomalous: for an order-of-magnitude reference — the exergetic $\eta_L \sim 3 \cdot 10^{-8}$ for the bacterium (§ 5.1), so the simulation $\eta_L \sim 10^{-5}$ lies in the range typical of informational efficiency. Drift without reset ($\tau_w \to \infty$, $T = 10^5$) exhibits an initial peak in the region $t \sim \tau_E$, after which $\eta_L$ decreases monotonically by 3 decades to $\sim 10^{-7}$ — substantively, the observer averaged $\hat P$ over $\sim 100$ independent regimes and almost does not predict the current $P(t)$. Drift with reset ($\tau_w \in \{100, 500, 2000\}$) gives a curve stationary on average without a long-term trend.

*Fig. 2: $\langle \eta_L \rangle$ as a function of $\tau_w$* (`paper/figs/Fig2.pdf`). A dome-shaped curve with a maximum $\langle \eta_L \rangle^* \approx 1.3 \cdot 10^{-5}$ at $\tau_w^* \approx 200$ steps; the position of the maximum $\tau_w^* \approx 0.2\,\tau_E$ agrees with (12) up to a logarithmic factor.

*Fig. 3: $\nu(t)$* (`paper/figs/Fig3.pdf`). In the drift-without-reset regime $\nu(t)$ grows monotonically from $\sim 0.3$ at $t \sim \tau_E$ to the asymptote $\nu_\infty \approx 0.999$ — almost complete nostalgia (the model accumulated over all past regimes is useless for predicting the current one); in the reset regimes $\nu(t)$ stabilises at the level $0.7$–$0.9$ depending on $\tau_w$. The empirical $\nu_c^{\text{emp}} \approx 0.44$ (at the half-life threshold of $\eta_L$) agrees with the structural statement of § 5.3 about the existence of $\nu_c \in (0, 1)$, but does not claim to be a universal quantitative estimate.

The full parameters of the three regimes and the $\tau_w$ grids, the pseudocode of the simulation, the KL-detuned proxy of the MI estimator with a consistency check on (Treves and Panzeri 1995; Kraskov et al. 2004), a discussion of the dome-shaped dependence as a bias–variance trade-off, the divergence $C_v^{\text{static}}/C_v^{\text{predictive}}$ as a diagnostic signal of nostalgia, the empirical $\nu_c^{\text{emp}}$ and its status — Supplementary § S6.

### 6.2 OU simulation: continuous drift of the logits

The OU simulation at $\sigma = 0.1$, $\lambda = 10^{-3}$ ($\sigma^2/\lambda = 10$, stationary variance $\sigma^2/(2\lambda) = 5$, $k = 8$, $K = k(k-1) = 56$, $T = 2 \cdot 10^4$ steps, 50 realisations) is a numerical illustration of Lemma 2 in the canonical slow-drift class of § 5.2: continuous OU drift of the logits of the transition matrix without resetting the posterior, in contrast to the PSP surrogate of § 6.1. The learner is a Robbins–Monro online approximation of the Bayesian update over the OU latent with a decreasing step $\eta_t = c_0/\sqrt t$, $c_0 = 3$ (asymptotically Bayes-optimal for regular parametric models Clarke and Barron 1990). The slow-driving condition $\lambda \tau_{\text{relax}} \ll 1$ is satisfied with a large margin; the OU-concentration condition ($\sigma^2/\lambda \ll 1$) is not satisfied here, so the conjecture of Supplementary § S8.1 is not tested in this parameterisation.

*Fig. 4: $\eta_L(t)$ for the OU regime* (`paper/figs/Fig4.pdf`). A monotonic fall of $\eta_L$ by 2–3 decades from $\sim 10^{-3}$ at $t \sim \tau_E = 10^3$ to $\sim 10^{-5}$ at $t = 2 \cdot 10^4$. Both budget variants ($N_{\max} = Rt$ and the full $N_{\max}^{\text{full}}$ with the summand $\nu(t)\cdot K/\tau_{1/2}$, $\tau_{1/2} = 100$) give qualitatively identical curves.

*Fig. 5: $\eta_L \cdot t$ vs $\ln(\lambda t)$* (`paper/figs/Fig5.pdf`). The numerical dependence exhibits a monotonic decay over the entire window $t \in [10^3, 1.5 \cdot 10^4]$; the late-time nostalgia $\nu(t > 10^4) = 0.82 \pm 0.03$ is high but does not reach the $\nu \to 1$ of the PSP surrogate. The numerical estimate $c \approx 0.02$ from $(*)$ gives the Lemma 2 prediction $\liminf \nu^{\text{theor}}(t) \ge 1 - c \approx 0.98$. (In this OU simulation the memory does *not* grow: $\dot E_{\text{grow}} = 0$, the number of logits $K$ is fixed, so $|M_{\le t}| = K$ — a property of the class of § 6.2, not a fitting of the denominator $(*)$ to the observed $\nu$.) The observed value $\nu(t > 10^4) = 0.82 < 1 - c$ lies *below* the asymptotic $\liminf$ floor (8a); this is consistent with the bound, not a violation of it: at finite $t$ the under-vanished stale fraction contributes a residual predictiveness, and finite-$t$ values lawfully lie below the asymptotic $1 - c$. The gap $\approx 0.16$ is explained by the finiteness of the observation time and the approximate fulfilment of per-bit-uniformity (Remark 5).

The OU simulation reproduces the collapse regime of § 5.1 on the canonical class of § 5.2 with continuous drift — without resetting the posterior, in contrast to PSP.

The full simulation parameters, the Robbins–Monro online learner and its asymptotic Bayes-optimality, the numerical estimate of $c$ for Robbins–Monro and the divergence from prediction (8a), the wide-window negative control of the adiabatic scan — Supplementary § S6. The full reproducibility parameters (seed `SEED = 20260524`, directories, library versions) — Supplementary § S4.

### 6.3 Conclusions from the simulations

The two simulations — PSP (§ 6.1) and OU (§ 6.2) — perform **four concrete functions** with respect to the formal theory of §§ 4–5. They do **not** claim empirical validation on real biological systems (this is the task of the pre-registered protocols of § 8.6) and do **not** test the conjecture of Supplementary § S8.1 on the adiabatic asymptotics of $L_{\text{excess}}$ (for which the OU-concentration condition $\sigma^2/\lambda \ll 1$ is needed, satisfied neither in PSP nor in OU § 6.2; the corresponding parametric scan — Supplementary § S8.3).

**(1) Structural confirmation of the three regimes of § 5.1.** PSP reproduces all three qualitatively distinct regimes of $\eta_L(t)$ — stationary (a horizontal asymptote $\sim 6 \cdot 10^{-5}$), collapse (a monotonic fall by 3 decades under drift without reset), and reset (a curve stationary on average at finite $\tau_w$). OU reproduces the collapse regime on the canonical class of § 5.2 with *continuous* drift, without resetting the posterior. The structural prediction of § 5.1 (a three-component classification by the relation of $\dot{M}_{\text{refresh}}$ and $\dot{M}_{\text{grow}}$) is *qualitatively confirmed* in both parameterisations.

**(2) Numerical check of Lemma 2 in its domain of applicability.** PSP with drift without reset gives $\nu_\infty \approx 0.999$ — almost complete nostalgia, consistent with (8a) in the regime $c \to 0$ (when updating has ceased relative to $\tau_E$). OU at $c \approx 0.02$ gives $\nu(t > 10^4) = 0.82 \pm 0.03$ against the predicted $\liminf \ge 1 - c \approx 0.98$; the gap $\approx 0.16$ is **below** the theoretical bound. This is **not a refutation** of Lemma 2 (it is a $\liminf$ at infinity), but a quantitative characterisation of the finite-$t$ correction and the approximate fulfilment of per-bit-uniformity (Remark 5 § S5.1) for the Robbins–Monro learner with a decreasing step $\eta_t = c_0/\sqrt t$. The wide-window negative control (Supplementary § S6) shows that under a "wrong" choice of the analysis window this gap can grow to $\sim 1$, which **strengthens** the methodological discipline of the protocol § 5.2 item 5 (the ban on post hoc redefinition of the object of measurement): the same numerical artefact may be read either as a refutation or as a confirmation depending on the analysis, which is why the fit window is pre-registered.

*An honest caveat about the quantitative gap.* Three observations are tied into one coherent picture. (i) The vanishing of the conditional MI of the stale fraction (§ 5.3, step (i)) is asymptotic: it sets in over *several* $\tau_E$, not in one step ($\rho = e^{-\Delta t/\tau_E}$ at $\Delta t = \tau_E$ leaves a residue of several bits, not zero), so over a finite observation window part of the formally stale bits still retains nonzero informativeness. (ii) The gap $\nu = 0.82$ vs $\liminf \ge 0.98$ ($\approx 0.16$) at $T = 20\,\tau_E$ is explained jointly by the finiteness of the observation time and the saturation of softmax in the regime of strong logit fluctuations ($\sigma^2/\lambda = 10$ in § 6.2), where the decay of correlations to the per-bit level is slowed relative to the idealisation (8a). (iii) Therefore the simulations confirm Lemma 2 *qualitatively* (monotonic $\nu \to 1$ in the limit $c \to 0$, existence of a positive asymptotic floor), whereas the quantitative check of the bound $1 - c$ as a point prediction is deferred to the pre-registered protocol on E. coli (§ 8.6) with controlled $\tau_E$ and a weak fluctuation regime.

**(3) Existence of $\nu_c$ and $\tau_{\text{reset}}^*$ — operationalisation of the structural results of §§ 5.3 and (12).** The dome-shaped curve $\langle\eta_L\rangle(\tau_w)$ (Fig. 2) with a maximum $\tau_w^* \approx 0.2\,\tau_E$ is an *operational instantiation* of (12) up to a logarithmic factor; the existence of an optimum $\tau_{\text{reset}}^*$ as a function of the structural parameters ($\tau_E$, $\dot{M}_{\text{refresh}}$, $|M_{\le t}|$) is not a postulate but a measurable effect on a synthetic benchmark. The empirical threshold $\nu_c^{\text{emp}} \approx 0.44$ (at the half-life of $\eta_L$) is a concrete numerical value consistent with the structural statement of § 5.3 about the existence of $\nu_c \in (0, 1)$; it does *not* claim to be a universal critical exponent (this is an open problem § 8.5) but gives an **operational handle** for designing the pre-registered prediction P.1 (§ 8.6) — the divergence $C_v^{\text{static}}/C_v^{\text{predictive}}$ as a diagnostic signal of the regime threshold.

**(4) Bridge to the ML operationalisation of Lemma 2.** The OU simulation with the Robbins–Monro learner is a structurally equivalent online continual-learning benchmark with a continuously drifting target distribution [discussion and a controlled paradigm case on the permuted-digits mini-analogue — § 8.4]. The numerical estimate $c \approx 0.02$ from $(*)$ for the specific parameterisation shows that **the operational computation of $c$ is possible** for systems with a standard memory-refresh scale; here $c$ is an order-of-magnitude estimate requiring empirical calibration. An important caveat: the naive ML form of the prediction (growth of $\nu^{\text{op}}$ with $\rho_p$) is *refuted* by the controlled paradigm case of § 8.4 as an artefact of underfitting — a small MLP under online SGD does not cleanly instantiate the law $c \to \nu$; the load-bearing falsifier of the carry-over of Lemma 2 remains biological (E. coli, § 8.6).

**What the simulations do *not* do.** (a) They do not validate the theory on real biological systems — this is the task of § 8.6 P.2 (E. coli chemotaxis). (b) They do not prove the conjecture of Supplementary § S8.1 about the asymptotics $(K/2)\ln(\lambda t)$ for $L_{\text{excess}}$ — for this the adiabatic OU regime $\sigma^2/\lambda \ll 1$ is needed, which is tested separately in Supplementary § S8.3 with a parametric scan (status: open conjecture, not a theorem). (c) They do not close the question of the carry-over of Lemma 2 to non-OU classes (PSP with Poisson reset, multi-scale drift, fractional OU) — Remark 4 § S5.1 explicitly restricts the rigorous proof to the OU class.

Summary: the simulations **confirm the qualitative structure** of the theory (three regimes, the existence of $\nu_c$ and $\tau_{\text{reset}}^*$, $\liminf \nu \ge 1-c$ as a lower bound), **quantitatively characterise** the finite-$t$ corrections and the conditions of applicability of per-bit-uniformity, **give an operational handle** for the pre-registered predictions (P.1, P.2 § 8.6; ML prediction § 8.4), and **remain within the domain of applicability** of the theory through explicit fulfilment of the conditions $(*)$ and slow-driving.

## 7. Connection with active inference and unsupervised learning

The connection between $\eta_L$ and Friston's free energy principle (Friston 2010; Friston et al. 2017a) is formulated in (Andriishin 2026, § 4.4): the requirement of self-payment turns the FEP into a discriminative criterion, and $\eta_L$ operationally closes the gap between the general formalism and the requirement that the energy belong to the modelling loop. In the non-stationary regime the Pearl/Friston blanket demarcation (Bruineberg et al. 2018; Bruineberg et al. 2022) is preserved: the correspondence between $\eta_L$ and EFE is meaningful only for a Friston blanket; for systems with externalised payment $\eta_L(t) \to 0$ — there is no correspondence.

In the non-stationary regime $\eta_L(t)$ structurally corresponds to the epistemic component of the expected free energy $G[\pi]$ (Friston et al. 2017a; Parr et al. 2022) in the slow-driving limit. **This correspondence is fixed as a *heuristic correspondence* (Supplementary § S7.1), not as a proved inequality**; the conjecture about an equivalent inequality (Supplementary § S7.2) is an open problem. Categorial caveat: the epistemic value $I_q(s'; o' \mid \pi)$ of AIF and $I(X_t; X_{t+\tau} \mid \hat P(t))$ of the present work are different MI objects (policy-contextual vs observational); the bridge through the sufficiency of $\hat P(t)$ gives a lower bound, not an equality (Supplementary § S7.3). The correspondence is *partial in the epistemic component*, not an equivalence; the pragmatic component through the preference $C$ remains outside the thermodynamic scale. A *hypothetical* active-inference algorithm with an extended objective function including $\dot E_{\text{store}}$ in the pragmatic part of the preferences would balance the reduction of nostalgia and the exploitation of the model with the optimum (12); existing AIF implementations (Friston et al. 2017a; Parr et al. 2022; Friston et al. 2021) optimise $G[\pi]$, not $\eta_L(t)$. The connection with intrinsic motivation (Schmidhuber 1991; Friston et al. 2017b) and meta-learning (Finn et al. 2017): in a non-stationary environment exploration is thermodynamically necessary for $\eta_L > 0$.

*Informational rejuvenation metric.* The ratio of the refresh and growth rates of memory

$$\rho(t) = \frac{\dot{M}_{\text{refresh}}(t)}{\dot{M}_{\text{grow}}(t)} \tag{13}$$

is a diagnostic independent of the direct measurement of $\nu(t)$: $\rho \gg 1$ — updating outpaces growth, nostalgia is low; $\rho \ll 1$ — the archival layer becomes nostalgic. The empirical calibration of $\rho_c$ is a separate study.

The full formal correspondence of $\eta_L$ and EFE with unit conversion through $\beta = 1/(k_B T)$, the heuristic correspondence for the increment $\Delta G^{\text{epist}}$ vs. $\Delta I_{\text{pred}}$, the conjecture about an equivalent inequality, the categorial caveat $I_q(s'; o') \ne I(X_t; X_{t+\tau})$, the context of sophisticated inference (Friston et al. 2021), the rebuttal of Andrews (2021) / Williams (2020) through a cross-reference to (Andriishin 2026, § S4.4a), and the detailed epistemic–pragmatic interpretation of $\rho(t)$ — Supplementary § S7.

## 8. Discussion, open problems, and conclusion

### 8.1 Summary

The work methodologically extends (Andriishin 2026) to the non-stationary regime: three structural results, a methodological toolkit, one open quantitative conjecture.

**Structural results.**

1. **Cumulative $N_{\text{max}}(t)$ and Lemma 1.** A cumulative budget with a decomposition into operational payment,
   holding, and memory expansion (§ 2); Lemma 1 on the minimal storage cost through cyclic rewriting
   in the naive-refresh regime with polynomial $1/\varepsilon$ and a logarithmic majority-vote alternative.

2. **Extension of predictive information and Lemma 2.** The one-step $I_{\text{pred}}(t,\tau) = I(X_t; X_{t+\tau}
   \mid \hat P(t))$ on $M_t \cup M_{\lt t}$ through the chain-rule decomposition (§ 3). Nostalgia with a two-level
   normalisation $\nu^{\text{theor}}$/$\nu^{\text{op}}$ (§ 4.1); Lemma 2: $\liminf_{t\to\infty}\nu^{\text{theor}}(t)
   \ge 1 - c$ at $c \lt 1$, with $c$ expressed through the refresh rate and $\tau_E$ (§ 4.4). $\eta_L(t)$ with three
   regimes (§ 5.1).

3. **The regime threshold $\nu_c$ and the optimal $\tau_{\text{reset}}^*$.** The postulate of a regime threshold
   "adaptive memory — nostalgic collapse" at $\nu_c^{\text{theor}} \in (0, 1)$ as a balance condition
   (§ 5.3); the optimal interval of periodic reset $\tau_{\text{reset}}^*$ (12).

**Methodological toolkit**: (i) $\nu^{\text{theor}}$ vs $\nu^{\text{op}}$ (§ 4.1); (ii) $C_v^{\text{static}}$ vs $C_v^{\text{predictive}}$ as a diagnostic of nostalgia (§ 4.3); (iii) the protocol against ad hoc rescues (§ 5.2, Supplementary § S5.2 — the ban on substituting the object of measurement); (iv) pre-registered protocols P.1, P.2 (§ 8.6).

**The quantitative part in Supplementary § S8.** A conjecture about the adiabatic asymptotics $\dot{\eta}_L^{\text{excess}} \sim -K\ln(\lambda t)/(2Rt^2)$ for $L_{\text{excess}}(t)$ in the OU class at $\sigma^2/\lambda \ll 1$; the parametric scan supports the functional form, the convergence to $K/2$ is not rigorously proved. The deferral to the supplementary is a consequence of the discipline of item 5 of the protocol.

*The strongest objection and the response to it.* Let us name directly the most incisive objection to the work: the central adiabatic asymptotics (Supplementary § S8.1) did not pass two independent numerical tests — neither PSP (§ 6.1) nor the OU simulation (§ 6.2) — since both regimes do not satisfy the OU-concentration condition $\sigma^2/\lambda \ll 1$, and therefore the asymptotics $B \to K/2$ remains unproved (§ 8.5). This is fixed as an *open conjecture*, not a result. The response: the load-bearing core of the work *does not depend* on this failure. Lemmas 1 and 2, the computable constant $c$ and the lower bound $\liminf\nu^{\text{theor}}\ge1-c$, the separation $\nu^{\text{theor}}/\nu^{\text{op}}$, the self-payment criterion, and the methodological protocol — structural results proved or formulated independently of the adiabatic asymptotics; the latter is a *secondary* quantitative superstructure over $L_{\text{excess}}$, voluntarily deferred to the supplementary precisely because it did not reach the status of a theorem. The refutation of the $K/2$ asymptotics touches none of the core statements.

The numerical illustrations (§ 6.1 PSP, § 6.2 OU) demonstrate Lemma 2 and the three regimes of $\eta_L(t)$ in the canonical class of § 5.2; $\eta_L(t)$ corresponds to the epistemic component of EFE (§ 7); the informational rejuvenation metric $\rho(t)$ is introduced (§ 7).

**Structural results and context-dependent specifications.** *Structural results* (general statements independent of the choice of parametric class): self-payment on the non-stationary budget (§ 2.6); the stationary lemma as a special case; $\eta_L(t) \in [0, 1]$ in the irreversible regime; the Landauer bound (1); Lemma 1; Lemma 2; the existence of $\nu_c^{\text{theor}}$ (with the working assumption about stability under a change of normalisation, § 4.1); the separation $\nu^{\text{theor}}$/$\nu^{\text{op}}$; the protocol § 5.2 / Supplementary § S5.2. *Context-dependent specifications* (dependent on the specific parametric class or operational choice): the drift class (OU/PSP); the numerical estimate of $\nu_c^{\text{op}}$; the choice of window $\tau$; the operationalisations of $C_v$; the open quantitative conjecture of Supplementary § S8.

### 8.2 What is NOT done in this work

The work is limited to theory. Beyond its scope are — (i) the empirical calibration of $\eta_L(t), \nu(t), \rho(t)$ on specific biospheres, cities, LLMs (requires the operationalisation of $E_{\text{store}}, E_{\text{grow}}$ for each class); (ii) a full theory of cognitive systems as an $\eta_L$ problem (psychological memory, forgetting, everyday nostalgia as phenomenal states of a subject — a categorially distinct referent from the informational nostalgia of the present work); (iii) the philosophical consequences of distinguishing the current and archival layers of memory (semantic vs episodic, the role of forgetting in identity).

### 8.3 LLM-corp as a cross-disciplinary illustration (a borderline case)

The class of large language models trained on a static corpus and operated in a drifting environment figures as a *cross-disciplinary illustration* of the applicability of the formalism to non-biological systems, not as an independent empirical application (the central prediction of the work is *E. coli* chemotaxis, § 8.6 P.2).

*Boundary and self-payment* (§ 2.6). In a strictly stateless LLM-as-agent the self-payment loop is structurally unclosed: $\dot E_{\text{store}}$ is paid for by an external data centre, and $\eta_L$ is undefined. By the escalating boundaries of paper #1 (Andriishin 2026, § 3.4 + Supplementary § S3.4), $\eta_L(t)$ is defined only starting from the $B_4$ boundary (the corporation as a legal entity): the current layer $M_t$ is the KV-cache plus the operational state of the corporation, the archival $M_{\lt t}$ is the weights plus the corporate corpus, the environment $X_E$ is the query statistics drifting in vocabulary (Lazaridou et al. 2021; Dhingra et al. 2022). This is an *honest demotion*: the formalism itself removes unclosed configurations from the domain of applicability rather than being stretched onto them.

*Two-layer prediction.* For the archival layer of an overfitted LLM without continual pre-training ($\dot{M}_{\text{refresh}}^{\lt t} \approx 0$, $c \approx 0$), Lemma 2 gives $\liminf \nu_{\lt t}^{\text{theor}}(t) \ge 1$ — the weights age as the language drifts, which qualitatively amounts to a monotonic growth of perplexity on texts after the cutoff — a **known observation** (Lazaridou et al. 2021; Dhingra et al. 2022), not an original prediction. For the current layer of long-context models, $c$ for the $M_t$ layer $\gg 1$, and Lemma 2 *is not applicable* to it (the sliding-window regime, an explicit self-limitation through $(*)$).

*Refresh channels and an open question.* In-context learning updates only $M_t$; continual pre-training and RAG with an updatable index (Lewis et al. 2020; Borgeaud et al. 2022; Khandelwal et al. 2020; Izacard et al. 2023) are a refresh channel for $M_{\lt t}$ (comparison of perplexity with/without updatable RAG is an operational prediction). Open question: does mesa-optimisation (Olsson et al. 2022; Garg et al. 2022; Bricken et al. 2023) induce an implicit update of the archival layer through ICL; if not — $\nu_{\lt t}^{\text{theor}} \to 1$ remains for all long-context LLMs without RAG.

### 8.4 Connection with continual learning

Catastrophic forgetting (McCloskey and Cohen 1989; Goodfellow et al. 2014; Kirkpatrick et al. 2017; Lopez-Paz and Ranzato 2017; Riemer et al. 2019; Parisi et al. 2019) has a natural interpretation within the framework of the theory; a parallel apparatus is the stochastic thermodynamics of learning (Goldt and Seifert 2017). Concept drift (Gama et al. 2014; Bifet and Gavald`a 2007) is unified through the bridge (8): drift detectors trigger at $\nu(t) > \nu_c$. Catastrophic forgetting is an attempt to reduce $\nu(t)$ through direct rewriting without allocating $E_{\text{grow}}$; EWC, replay buffers, progressive networks — investment in $E_{\text{grow}}$.

Nostalgia and catastrophic forgetting are two poles of the balance of $\dot{E}_{\text{store}}$ and $\dot{E}_{\text{grow}}$: nostalgia is an *excess* of non-predictive memory, forgetting a *deficit* of prematurely erased memory. The optimum is (12). At $\lambda \to 0$ nostalgia reduces to *overfitting*.

*Contrast with renaming — comparison with ML precedents.* To the objection "nostalgia is a renaming of concept drift + forgetting + dynamic regret": the theory contains constructions structurally close to known ML results but operationally and scale-wise distinct from them.

*The protocol acted on the authors themselves.* The clause of the protocol (§ 5.2, Supplementary § S5.2) forbidding the post hoc redefinition of the object of measurement is a substantive contribution, not a methodological ritual: it acted precisely on the authors of the present work. After two numerical tests (PSP § 6.1 and OU § 6.2) failed to confirm the adiabatic asymptotics, a tempting "rescue" would have been to redefine the numerator — to replace the operationally measurable $I_{\text{pred}}(t,\tau)$ with the counterfactual $L_{\text{excess}}$, under which the asymptotics holds — while keeping the name of the theory. The protocol forbade this move: $L_{\text{excess}}$ is deferred to § S8 as a secondary metric, and the asymptotics is marked as an open conjecture. Without the protocol the same material could have been presented as a "confirmation". A discipline that works against the author's own theory is not rhetoric but a real constraint.

**(a) Lemma 2 vs Besbes–Gur–Zeevi 2015 (dynamic regret).** (Besbes et al. 2015) (Theorem 2) give, for a sliding-window learner, a regret $R_T = O(V_T^{1/3} T^{2/3})$; the parallelism is obvious. The difference: this is a *sample-level upper bound* on the regret at $V_T = o(T)$, whereas Lemma 2 is a *bit-level lower bound* on nostalgia with an energetic budget $N_{\max}(t)$ in the OU class of § 5.2, that is, a thermodynamic lower bound on a measure that sample-level regret does not reflect.

**(b) $\tau_{\text{reset}}^*$ vs optimal sliding-window length (Hazan and Seshadhri 2009; Karnin and Anava 2016).** There $\tau_{\text{reset}}$ is an *algorithmic* restart parameter without a thermodynamic interpretation; in (12) it is a *physical* interval of periodic memory wipe, optimising $\eta_L$ through the trade-off $\dot E_{\text{store}}$ vs $\dot E_{\text{grow}}$ (the biological instantiation is the generation time of epigenetic reprogramming, § 4.4).

**(c) $\nu_c$ vs ML phase transitions.** $\nu_c \in (0,1)$ (11b) is *not* equivalent to grokking (Power et al. 2022) (reorganisation of weights at a fixed environment) or double descent (Belkin et al. 2019) (generalisation vs capacity at a fixed task): $\nu_c$ is a transition in *memory composition* under a drifting environment and a fixed architecture; the connection is non-trivial (§ 8.5).

**(d) The thermodynamic scale $k_B T \ln 2$ is absent from the ML literature.** Nostalgia as an energetic quantity is not defined there; $\dot{E}_{\text{store}}^{\text{nostalg}} \ge k_B T \ln 2 / \tau_{1/2}$ is a bridge between the informational measure $\nu$ and the physical budget, unraisable within ML formalisms (Still et al. 2012; Gama et al. 2014; Kirkpatrick et al. 2017). In the biological paradigm case, by contrast, the energetics of adaptation is established: the price of the accuracy and rate of sensory adaptation is quantitatively described for *E. coli* chemotaxis (Lan et al. 2012; Mehta and Schwab 2012). The contribution of the present work here is not the first introduction of energetics into chemotaxis, but the *non-stationary lower bound on the accumulation of nostalgia* — on the flow of holding the ageing fraction of memory over time, not on a single adaptation.

*An operational handle on $c$ for ML — an illustrative estimate on the permuted-digits mini-analogue.* The permuted-digits protocol (prototype — Permuted MNIST (Goodfellow et al. 2014); in the reproducible paradigm case below — `sklearn.load_digits`, $8\times8$) with permutation frequency $\rho_p$ gives $\tau_E = 1/\rho_p$; for SGD with learning rate $\eta$ on $N$ parameters the effective memory-refresh rate is $\dot M_{\text{refresh}} \approx \eta N$ effective scalar updates per SGD step (treated order-of-magnitude as bits; SGD updates $\sim \eta N$ such quantities on *each* step, not per permutation-switch event), with $|M_{\le t}| \approx N$. Substitution into the canonical $(*)$ ($c := \dot M_{\text{refresh}}\,\tau_E/|M_{\le t}|$) gives
$$c \;\approx\; \frac{\eta N \cdot (1/\rho_p)}{N} \;=\; \frac{\eta}{\rho_p}$$
($c \propto 1/\rho_p$; the exact value is calibrated empirically — the expression $\eta N$ here serves only as an order-of-magnitude factor of the refresh rate, not a dimensionless fraction of memory; order-of-magnitude at $N \sim 10^5$, $\eta \sim 10^{-3}$). A *fast* switch (large $\rho_p = 10^{-2}$/step) gives $c = \eta/\rho_p = 0.1 \ll 1$ — memory lags behind the drift, $\nu \to 1 - c \approx 0.9$ (high nostalgia); $\rho_p = 10^{-3}$/step — the boundary $c \approx 1$; a *slow* switch (small $\rho_p = 10^{-6}$/step) gives $c = 10^3 \gg 1$ — refresh dominates, the sliding-window regime, Lemma 2 is *inapplicable* through $(*)$, $\nu \to 0$ (low nostalgia). The illustrative point of formula $(c)$ takes $\eta \sim 10^{-3}$, whereas the reproducible paradigm case below operates at $\eta = 0.1$, $c \gg 1$ — these are different parametric points: the paradigm case ($\eta = 0.1$, $c \gg 1$) illustrates precisely that SGD sits in the sliding-window regime. Naive reading: nostalgia grows with the switch frequency $\rho_p$ (fast drift — memory does not have time to update — high $\nu$), formally agreeing with § 8.3 ($c \to 0 \Rightarrow \nu \to 1$). Let us note the classical protocol at once: in the setting of (Goodfellow et al. 2014) each permutation is trained to convergence, and only then is a switch made — this is a *rare* effective task switch, that is, small $\rho_p$ and large $c$ (and not large $\rho_p$, as one might erroneously conclude). Moreover, classical catastrophic forgetting (Goodfellow et al. 2014) measures the *drop of accuracy on past tasks*, not the ablatability of units on the current task — this is not the same quantity as $\nu^{\text{op}}$, and they cannot be identified.

*The ML form of the prediction and its proxying character.* The naive ML form of the prediction is as follows: for a trainable MLP without EWC/replay the nostalgic fraction should grow with $\rho_p$ ($1 - \nu \sim c = \eta/\rho_p$, that is, $\nu^{\text{op}}$ grows monotonically with the permutation-switch frequency). Here $\nu^{\text{op}}$ is operationalised through the ablation fraction — the fraction of hidden units whose zeroing does *not* raise the current-task loss above threshold. We stress: the ablation fraction in ML is only a *proxy* of the operational $\nu^{\text{op}}$ of § 4.1 (which is canonically $1 - I_{\text{pred}}/C_\mu^{\text{emp}}$), not an identity; in ML it additionally requires control of learnedness (see below).

*Reproducible paradigm case and an honest negative result (`simulations/permuted_mnist/`).* The naive ML form was checked by a controlled paradigm case on the permuted-digits mini-analogue (`sklearn.load_digits`, $8\times8$, MLP $64\to32\to10$, online SGD, $\eta = 0.1$, $H = 32$, 4 seeds) — and **does not survive**. Under the naive ablation metric $\nu^{\text{op}}$ a monotonic growth with $\rho_p$ was indeed observed, but the controlled paradigm case showed that this is an **artefact of underfitting**: $\nu^{\text{op}} \approx 1 - \text{learnedness}$ (under fast drift the network does not have time to learn the task, zeroing any unit is trivially harmless, $\nu^{\text{op}} \to 1$). After an accuracy gate (only regimes where the task is learned, accuracy $\gg$ chance, are interpreted) and the definition of *true* nostalgia (a unit is useless now **and** was useful on a recent past permutation — which separates stale-retained units from dead/never-useful ones) the predicted effect **disappears and reverses**: in the learned regime $\nu^{\text{op}}_{\text{true}}$ *anti*-correlates with $\rho_p$ (Spearman $= -1.0$; on $N = 4$ points only the *sign* is substantive — the reversal of the naive sign — not the magnitude of the correlation). Over the four learned points ($\rho_p \in \{10^{-3}, 5\cdot10^{-4}, 3\cdot10^{-4}, 10^{-4}\}$, accuracy $\{0.88,\,0.91,\,0.955,\,0.989\}$): $\nu^{\text{op}}_{\text{true}} = \{0.025,\,0.035,\,0.141,\,0.402\}$, dead fraction $\le 0.09$; the point $\rho_p = 3\cdot10^{-3}$ (accuracy $0.78$) is excluded as underfitted. Physical reason: SGD is a high-refresh learner, constantly sitting in the large-$c$ regime (sliding window), so a small MLP does **not** cleanly instantiate the law $c \to \nu$; the dominant driver of true nostalgia here is re-specialisation under long dwell ($1/\rho_p$), not the lag of memory behind the drift frequency.

Substantively this is the **third honest negative test** of the article — after the check of the PSP bound (§ 6.1) and the OU class for the adiabatic asymptotics (§ 6.2). The methodological protocol § 5.2 (the priority of falsifiability over confirmation) acted here too: without control of learnedness, a confound-driven "confirmation" would have been retained as valid; the control prevented this. The load-bearing falsifier of Lemma 2 remains **biological** — bacterial *E. coli* chemotaxis (§ 8.6); the naive ML form is refuted by its own paradigm case and is not presented as a confirmation.

*Class-restriction and no-free-lunch.* The Wolpert no-free-lunch (Wolpert 1996) gives a fundamental justification for narrowing the class (§ 5.2): a universal theory of non-stationary learnability without an a priori restriction of the class is impossible; the slow-driving class is a necessary narrowing, not a methodological choice.

### 8.5 Open problems

*The conjecture about the adiabatic asymptotics (Supplementary § S8).* The main open theoretical problem is the rigorous proof or refutation of the conjecture of Supplementary § S8.1 about the asymptotics $\dot\eta_L^{\text{excess}} \sim -K\ln(\lambda t)/(2Rt^2)$ for the cumulative excess-loss $L_{\text{excess}}(t)$ in the OU class at $\sigma^2/\lambda \ll 1$. The numerical parametric scan of Supplementary § S8.3 supports the functional form ($R^2 \ge 0.9995$ in the fit $A\cdot t + B\cdot \ln(\lambda t) + C$) but does not prove the convergence of the coefficient $B \to K/2$: the monotonic growth of $B/(K/2)$ from 0.54 to 0.85 as $\sigma^2/\lambda$ goes from $10^{-1}$ to $10^{-3}$ is consistent with a hypothetical finite-σ correction $B(\epsilon) = (K/2)(1 - g(\epsilon))$, but a pre-registered confirmation protocol requires at least 5–6 scan points with an explicit analytical form of $g(\epsilon)$.

*Extension of Lemma 2 to rapidly growing memory.* The current Lemma 2 relies on the boundedness of $|M_{\le t}|$ through the constant $c$ from $(*)$. At $|M_{\le t}| \propto t^\alpha$ ($\alpha \ge 1$) the condition $c \lt 1$ becomes trivially satisfiable, but the qualitative picture of nostalgia changes: under sufficiently rapid growth of capacity the nostalgic fraction may remain a finite fraction. A formal extension of Lemma 2 to this regime is an open problem of non-stationary Bayesian models.

*Connection with fluctuation theorems.* The theorems of Crooks (Crooks 1999) and Jarzynski (Jarzynski 1997); the generalisation with feedback control (Sagawa and Ueda 2010) gives a bridge to the decomposition (1) through the measurement channel. An application to the balance of $E_{\text{store}}, E_{\text{grow}}, E_{\text{actual}}^{\text{curr}}$ may give finer inequalities than the Landauer bound.

*Operationalisation of $E_{\text{store}}$ for cultural and biospheric memory.* $\tau_{1/2}$ is quantitatively different (epigenetics, genome, synapses); in cultural systems the carrier is ambiguous (book, digital archive, neural network) — a separate study.

*Refinement of the regime threshold $\nu_c$.* The existence of a rigorous transition with universal critical exponents and a numerical estimate for the paradigm classes is an open question. To this also belongs the *stability of $\nu_c$ under a change of normalisation* (§ 4.1): a rigorous theorem on the monotonicity of the regime threshold under the DPI inequality $\nu^{\text{op}} \le \nu^{\text{theor}}$, justifying the transition from $\nu_c^{\text{theor}}$ to $\nu_c^{\text{op}}$ in the pre-registered protocol § 8.6 (P. 1).

*Regularity of the softmax parameterisation.* A full proof of the transfer of OU decorrelation onto the nonlinear softmax parameterisation of the transition matrix (Remark 6 § 4.4) is a technical extension through Lipschitz continuity in total variation on the part of the OU trajectory bounded with probability $1-\delta$; it gives the convergence $\nu^{\text{theor}}(t) \to 1 - c$ in the specific OU parameterisation of § 6.2 / Supplementary § S6.3 beyond the illustrative $\liminf$ inequality (8a).

### 8.6 Empirical predictions

The theory formulates two core empirical predictions (P. 1, P. 2).

*Prediction 1: existence of a regime threshold $\nu_c^{\text{op}} \in (0, 1)$.* In adaptive non-stationary systems with a finite memory-refresh rate, a critical level of nostalgia $\nu_c^{\text{theor}} \in (0, 1)$ should exist (§ 5.3); it is operationally diagnosed through the empirical threshold $\nu_c^{\text{op}}$ by the divergence $C_v^{\text{static}}/C_v^{\text{predictive}}$ under observed environmental drift (§ 4.3): at $C_v^{\text{static}}/C_v^{\text{predictive}} \to \infty$ the system has crossed $\nu_c^{\text{op}}$ (see the working assumption § 4.1 about the connection of $\nu_c^{\text{theor}}$ and $\nu_c^{\text{op}}$). Refutation: the absence of divergence of the two proxies in systems with deliberately impaired adaptation (for example, frozen LLM weights in a drifting query stream; fully methylated bacteria without active demethylase correction).

*Prediction 2: slowdown of the adaptation rate upon blocking the receptor-methylation refresh loop (E. coli chemotaxis, the core qualitative prediction).* The biological paradigm case of Lemma 2 in a bacterial system is the chemotaxis adaptation loop: the methyltransferase CheR covalently methylates the methyl-accepting chemotaxis receptors (MCPs: Tar, Tsr, Tap, Trg) on conserved glutamate residues; the methylesterase CheB actively demethylates them in the phosphorylated mode, closing the negative-feedback loop of precise adaptation (Yi et al. 2000; Alon et al. 1999; Sourjik and Berg 2002; Endres and Wingreen 2006). The methylated states of the receptors realise the $M_{\le t}$ of Lemma 2 (an internal state storing the statistics of the background ligand concentration); CheR/CheB-mediated re-methylation is $\dot M_{\text{refresh}}$.

*Protocol.* Suppression of the refresh loop is realised through a **hypomorphic CheR mutant with reduced activity** ($cheR^{\downarrow}$) or controlled reduction of *cheR* expression from a titratable (IPTG-inducible) promoter (Alon et al. 1999, Fig. 2) — this preserves precise adaptation as a function (unlike a $\Delta cheB$ knockout, which gives an immediate locked-tumbling phenotype through hyper-phosphorylation of CheY, unrelated to the gradual nostalgia accumulation of Lemma 2). The environment is a **classical chemotactic gradient**: periodic switching of the concentration of aspartate / α-methyl-DL-aspartate (MeAsp) for the Tar receptor, or serine for Tsr, in a flow microfluidic chamber with a known profile; constant concentration is the control environment. The characteristic environmental drift time $\tau_E = 10$–$60$ minutes (the biologically relevant scale of methylation refresh is the minute scale of adaptation Sourjik and Berg 2002; Endres and Wingreen 2006); the horizon of the experiment is 5–10 generations (the typical generation time of *E. coli* in M9 minimal medium is 40–60 min), which gives 50–300 drift cycles over the horizon.

*Primary endpoint.* The adaptation rate $r_{\text{adapt}}(t)$ — the inverse of the recovery time of the tumbling frequency or CCW bias to baseline after a step change in attractant concentration, measured by FRET on the CheY–P / CheZ pair or by single-cell tracking. By Lemma 2, at $\dot M_{\text{refresh}}^{cheR^{\downarrow}} \ll \dot M_{\text{refresh}}^{WT}$ the nostalgic fraction accumulates under periodic concentration switching, and $r_{\text{adapt}}^{cheR^{\downarrow}}(t)$ slows relative to $r_{\text{adapt}}^{WT}(t)$ after several $\tau_E$ cycles. The *effect size* is formulated as a **comparison of the adaptation rate under periodic concentration switching vs constant concentration for the same pair of $cheR^{\downarrow}$ vs WT populations**: a genotype × environment interaction is predicted ($\Delta r_{\text{adapt}}^{cheR^{\downarrow}} [\text{periodic} - \text{constant}] / \Delta r_{\text{adapt}}^{WT}[\text{periodic} - \text{constant}] \lesssim 0.5$ after 5–10 switch cycles). Refutation: the absence of an interaction (ratio $\gtrsim 0.8$) given a passed test of membership in the class of slow drift of the attractant concentration.

*Alternative constructs and confounds.* (i) The $\Delta cheB$ knockout is excluded as a clean test of Lemma 2: constant immediate hypermethylation gives a locked-tumbling phenotype unrelated to gradual nostalgia accumulation. (ii) The knockout of the DNA adenine methyltransferase *dam* ($\Delta dam$) is an irrelevant confound: a mutator phenotype through impaired mismatch repair, not reducible to nostalgia. (iii) Switching the carbon source (glucose ↔ glycerol) is not a paradigm-case protocol: these substrates are not catabolites of Tar/Tsr and do not directly induce the CheR/CheB adaptation loop.

**The protocol against ad hoc rescues for P. 1 and P. 2** is parallel to § 5.2 / Supplementary § S5.2:

1. A pre-registered protocol (OSF/Zenodo, timestamp ≤ $T_0$) with a power analysis for the chosen test: for P. 1 —
   statistical discrimination of the two proxies with a fixed threshold; for P. 2 — comparison of the monotonicity
   of $r_{\text{adapt}}(t)$ between groups through a test of trends or an equivalent.
2. A pre-registered test of membership in the class of slow drift of the source (CUSUM and Ljung–Box, $p_{\text{class}}
   \ge 0.1$). At $p_{\text{class}} \in [0.01; 0.1]$ — a "grey zone": no more than one pre-registered repeat
   with an increased $N$ (Holm–Šidák); after two "grey" outcomes — item 4.
3. A ban on ad hoc appeals ("wrong taxon / wrong pharmacology / unrepresentative substrate") after a failure;
   mitigating factors are fixed before data collection.
4. A change of the theory's name upon systematic rejection: $\nu$-theory v2 with retention of the v1 predictions as
   a boundary condition. Retention of the name with a change of class is an inadmissible rescue through a change of the object of measurement.

Refutation is a qualitatively opposite observation under the satisfied conditions.

**Excess predictions relative to (Andriishin 2026).** Four classes of predictions not derivable from the stationary formulation:

1. Lemma 2 on the inevitability of nostalgia with an explicit computable lower bound $\liminf \nu^{\text{theor}}(t) \ge 1-c$
   through $c$ (§ 4.4).
2. The critical $\nu_c^{\text{theor}} \in (0, 1)$ separating the adaptive regime and collapse; the operational
   $\nu_c^{\text{op}}$ — a numerical estimate § 8.5, the working assumption about stability § 4.1.
3. The optimal $\tau_{\text{reset}}^*$ (12), correlating with the generation time of biological taxa.
4. The distinction between $C_v^{\text{static}}$ and $C_v^{\text{predictive}}$ as a diagnostic signal of nostalgia.

Additionally: the open conjecture about the adiabatic asymptotics $\dot\eta_L^{\text{excess}}$ (Supplementary § S8). The refutation of any of the four core statements refutes the non-stationary extension, not the stationary theory. The theory makes four excess empirical predictions relative to (Andriishin 2026); their independent verification sets the criterion of the progressiveness of the extension.

### 8.7 Conclusion

**Main conclusion:** *memory is not a record but an active thermodynamic process that continuously pays for itself; for any finite refresh rate and slow OU drift of the environment, part of the memory inevitably ages and keeps requiring energy without helping to predict* — this effect is quantified by Lemma 2 ($\liminf \nu^{\text{theor}} \ge 1 - c$ in the OU class) and is pre-registered for testing on bacterial chemotaxis (§ 8.6 P.2).

Real learning systems — bacteria, biospheres, neural networks, corporations — operate in a non-stationary environment with growing memory; the extension of $\eta_L$ (Andriishin 2026) requires accounting for the cost of storage and memory expansion, the concept of nostalgia, and the recognition of its inevitability in the OU class of slow drift. $\eta_L(t)$ decreases under drift; sustaining $\eta_L > 0$ requires active memory reset as a thermodynamic investment. Every held bit carries a payment by Lemma 1 (for feedback memory — in the form (2')); a nostalgic bit is a payment without gain.

Self-payment as a condition of applicability of (1) is an operational reformulation of the autonomy criterion in contemporary philosophy of biology: the capacity of a system to pay itself for holding its own model of the environment is quantified by the measurable count $-dF_X/dt$ over the window $\tau_d$, without imposing ontological commitments about the form of causal closure (§ 2.6, § 4.7). This operationalisation is consistent with Cleland's 2019 argument for a transition from lexicographic criteria of life to a *theory* of life (Cleland 2019) resting on testable physical invariants, and is in principle transferable beyond terrestrial biochemistry — the connection with the astrobiological search for biosignatures is developed in (Andriishin 2026, § 6) for the stationary regime; the non-stationary extension (periodic memory reset as a signature of autonomy on aeonic windows) is an open problem § 8.5.

---

## Appendix A. Normalisations and units of measurement

Throughout the work, $I_{\text{pred}}(t,\tau)$, $\nu(t)$, $\eta_L(t)$, $N_{\text{max}}(t)$, and related quantities are measured in **bits** (consistent with (1), where the denominator contains $k_B T \ln 2$). This convention reconciles the dimensions with the bit form of the Still–Sivak–Bell–Crooks inequality (Still et al. 2012, eq. 3–6) and with the usual scale of the Landauer bound. For conversion to nats, the factor $\ln 2$ is used: $I_{\text{[nat]}} = I_{\text{[bit]}} \cdot \ln 2$.

In Supplementary § S8 the cumulative excess-loss $L_{\text{excess}}(t)$ — a theoretical secondary metric, kindred to the information regret of a Bayesian learner in the sense of (Rissanen 1986; Clarke and Barron 1990) — is measured in **nats** (as in the original BNT literature (Bialek et al. 2001, § VI) and the MDL formulation); the explicit conversion $L_{\text{excess}}^{\text{[bit]}}(t) = L_{\text{excess}}^{\text{[nat]}}(t)/\ln 2$ is applied when comparing with the main text.

---

## Author Contributions

Conceptualisation, methodology, formal analysis, writing — preparation of the original draft, writing —
review and editing: A.A. The author has read and agreed to the published version of the manuscript.

## Funding

The research received no external funding.

## Ethics declarations

Ethical approval, consent to participate, and consent to publish are not applicable: the present work is a
theoretical study with numerical simulations that did not involve humans, animals, biological material, or
related data as subjects.

## Data availability

No primary empirical data were generated for this work. The reproducible code of the numerical simulations (Markov chains with drift § 6.1; the controlled ML paradigm case on the permuted-digits mini-analogue § 8.4 in `simulations/permuted_mnist/`) is openly available in the GitHub repository `andriishin/landauer-nostalgia-oa` (https://github.com/andriishin/landauer-nostalgia-oa) and is archived on Zenodo together with a preprint of this manuscript (DOI: 10.5281/zenodo.20653051).

## Use of AI Tools

In accordance with the Springer Nature editorial policy on AI tools in scholarly writing, AI is not listed as a co-author and the author retains full responsibility for all content of the manuscript. In preparing the present manuscript and supporting materials, the author used Anthropic Claude (Opus 4.6, 4.7, 4.8) for five tasks. (i) Generation and verification of reproducible Python code for the numerical simulations in `simulations/` (§ 6.1, § 6.2, § S8). (ii) Cross-checking bibliographic metadata against local copies of the cited sources. (iii) Editorial revision of human-generated text at the level of readability, grammar, and style (copy-editing); per the Springer Nature guidance, such use does not require a separate declaration but is stated here for full transparency. (iv) Translation of the manuscript and supporting materials (including the READMEs in `simulations/`) from the Russian primary source into English — the journal's submission language — with subsequent terminological and stylistic revision; the author checked the translation and bears full responsibility for it. (v) Assistance with the derivation and verification of the formal apparatus (the proofs of Lemmas 1–2, the OU estimates, and the explicit constants) under the author's direction; the author independently checked all proofs and bears full responsibility for them. All scientific statements — the formulation of the non-stationary $\eta_L(t)$, Lemma 1 on the minimal storage cost, Lemma 2 on the inevitability of nostalgia in the slow-drift class, the methodological refinement of $C_v$, the open conjecture about the adiabatic asymptotics $\dot{\eta}_L$ — were conceived by the author; their formal write-up was carried out with AI assistance under the author's control. The AI tools did not autonomously generate scientific content. The author checked every AI-assisted fragment and accepts full responsibility for the content of the manuscript and supporting materials.

## Competing interests

The author declares the absence of any relevant financial or non-financial interests that could be
perceived as influencing the content of this article (over the three years prior to manuscript submission).

## References

The reference list is generated automatically from `paper/refs.bib` at LaTeX compile time using a
Springer **Author-Year** bibliography processor (`sn-mathphys-ay.bst` or the journal-specific
Author-Year style — verify at submission). In-text citations follow the name-and-year convention
(`Author Year`; `Author and Coauthor Year`; `Author et al. Year`), and the reference list is ordered
**alphabetically by the first author's surname**, as required by *Theory in Biosciences* submission
guidelines (`docs/biosciences-requirements.md` § 8).
