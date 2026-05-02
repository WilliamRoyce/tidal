# Speaker notes — practice talk (2026-04)

Raw research material for the spoken script. One heading per slide. Phase D turns these into `\note{}` blocks in the Beamer source; Phase E turns them into continuous prose.

The talk now follows the project abstract's spine: standard Gertsenshtein is astrophysically useless → systematic phenomenological scan over *all parity-even quadratic curvature + torsion terms* → derive PDEs symbolically → screen for causality-violating regions via all-mode eigenvalue stability before sampling → identify any subsector that can amplify.

**Honest framing note for the causality slide**: the abstract pitched "characteristic speeds for causality checking" as a method, but in practice the project never implemented characteristic-speed analysis in the formal PDE sense. What it *did* implement, and what fills the same role, is eigenvalue stability pre-screening: compute the largest real part of any eigenvalue across all spatial modes, reject parameter points with growing modes. This is the working causality gate. Mention the abstract's framing only to say what the implementation evolved into.

---

## Narrative spine (one sentence)

> "The standard Gertsenshtein effect is too weak to ever observe; I'm scanning the space of parity-even quadratic extensions of GR — curvature-squared, torsion-squared, and their cross-couplings — to find any subsector that amplifies it, and to check whether the surviving theories are causal."

Every slide either sets up or executes a clause of this sentence.

---

## 1. Title
- One sentence intro: "I'm a master's student with the Handley group, working on whether modifying gravity can amplify a particular old graviton-photon conversion process. Practice talk — feedback welcome."

## 2. Motivation — *why anyone should care*

**The headline numbers**:
- For a magnetar field $B \sim 10^{15}$ G over its $\sim 10$ km radius, the Gertsenshtein conversion probability is $P \sim 10^{-10}$ — completely negligible against any astrophysical photon background.
- For a laboratory-scale $B \sim 10$ T, $D \sim 1$ m, $P \sim 10^{-40}$ — there is no foreseeable detector that can see this.
- **The standard Gertsenshtein effect is astrophysically useless and laboratory-impossible.** That isn't an opinion; it's the textbook answer for Einstein-Maxwell.

**So the question is forced**: *if Einstein-Maxwell can't give us a useful signal, what can?* New physics is required for any significant amplification.

**My approach** (paraphrasing the abstract): unapologetically phenomenological. Don't pick one theory and defend it — instead, systematically scan the space of *parity-even quadratic curvature and torsion terms* added to the Einstein-Maxwell Lagrangian, and ask which (if any) amplify the conversion.

**Three-beat promise** for the rest of the talk:
1. The channel — what is the Gertsenshtein effect, why is it weak, why might modifying gravity help?
2. The survey program — what's the parameter space we're scanning, and how do we scan it?
3. What we've learned so far — partial map of the survey, what amplifies, what doesn't.

---

## 3. The Gertsenshtein effect

**Setup to say out loud:**
- Einstein-Maxwell linearised around flat Minkowski + uniform background $B_0$ along $x$ (Coulomb gauge $\bar A_y = -B_0 z$).
- Expand metric perturbation $h$ (graviton) and EM perturbation $a$ (photon) to first order.
- TT gauge on $h$, Lorenz on $a$, plane-wave reduction.

**Coupled linearised system**:
$$\partial_t^2 h_+ = \partial_z^2 h_+ - \tfrac12\kappa^2 B_0^2 h_+ - \kappa^2 B_0\,\partial_z a_y$$
$$\partial_t^2 a_y = \partial_z^2 a_y + B_0\,\partial_z h_+$$

- Asymmetric coupling: gravity kinetic normalisation $1/(2\kappa^2)$ vs Maxwell $1$.
- Beat frequency $\Delta\omega = \kappa B_0$ — geometric mean of the two asymmetric couplings.

**Headline formula**: $P(g\to\gamma) = \sin^2(\tfrac12\kappa B_0 D)$.

**Physics intuition** (the analogy): same math as two-state QM mixing, same as axion-photon Primakoff. The deeper picture: the background $B_0$ contributes to the *graviton's* effective dispersion via its stress-energy; the photon stays massless. The detuning sets the beat frequency.

---

## 4. Why standard Gertsenshtein is useless — and what the literature says now

**Why it's so weak**: $\kappa \sim 10^{-19}$ GeV$^{-1}$ — the Planck-scale suppression of gravitational interactions sits inside the conversion frequency. Doubling $B_0$ doubles the conversion rate; you need 20 orders of magnitude.

**Literature timeline**:
- **Gertsenshtein 1962** — original paper, gave the formula, recognised it was hopelessly small.
- **Boccaletti et al. 1970** — rederived with plasma corrections.
- **Raffelt & Stodolsky 1988** — unified framework with axion-photon mixing.
- **2020s revival** — high-frequency GW detection programme: Aggarwal et al. *Living Reviews* 2021, Berlin et al. 2022 (microwave cavity haloscopes), Domcke & Garcia-Cely 2021 (radio-telescope HFGW detection).

**Why the revival now**: HFGW detector schemes use *exactly* the inverse Gertsenshtein process to convert hypothetical HFGW signals into photons in a cavity. So the small standard-GR conversion rate sets the noise floor for the entire HFGW programme — and any BSM physics that *amplifies* the conversion is *directly* relevant for detection.

**Worth dropping** (`gertsenshtein_formula.tex:26-95`): the formula is *still* being rederived in the 2020s, and modern papers disagree on prefactors of $\sqrt{4\pi}$ depending on graviton normalisation conventions (Palessandro & Rothman 2023 vs Dandoy-Lella et al. 2024). That's a fingerprint of an under-explored corner of theoretical physics.

---

## 5. Beyond GR: Poincaré gauging

**The framing the audience knows**: GR uses a torsion-free, metric-compatible connection (Levi-Civita). That's a *choice*, not a consequence.

**Utiyama 1956, Kibble 1961, Sciama 1962**: gravity as a gauge theory of the Poincaré group, analogous to gauging $U(1)$ for QED. Two subgroups, two field strengths:

| Gauge subgroup | Connection         | Field strength                    |
|----------------|--------------------|-----------------------------------|
| Lorentz        | spin connection $\omega$ | **curvature** $R^{ab}{}_{\mu\nu}$  |
| Translations   | tetrad $e$         | **torsion** $T^a{}_{\mu\nu}$        |

**Tagline**: *"GR keeps one of the two field strengths of the Poincaré group. PGT keeps both."*

**Geometric intuition** (loop cartoon — placeholder image): parallel-transport around a closed loop has *two* mismatches — rotational (curvature) and translational (loop fails to close — torsion).

**Why take it seriously**:
- Couples naturally to fermion spin (Einstein-Cartan, Hehl et al. 1976).
- Minimal extension — smallest well-motivated deviation from GR.
- Gives the most parameters to play with for *amplifying physics*, which is what we want.

---

## 6. The systematic landscape — *all parity-even quadratic terms*

**The phenomenological scan** (this is the slide that does the most work for the abstract):

The Lagrangian we consider is
$$\mathcal{L} = \frac{R}{2\kappa^2} + \alpha_1 I_1 + \alpha_2 I_2 + \alpha_3 I_3 + \beta_1 R^2 + \beta_2 R_{\mu\nu}R^{\mu\nu} + \beta_3 R_{[\mu\nu]}R^{[\mu\nu]} + \dots + \delta_1 R_{[\mu\nu]}F^{\mu\nu} + \dots - \tfrac14 F_{\mu\nu}F^{\mu\nu}$$

The terms split into three families:

1. **Torsion-squared invariants** — three of them ($I_1$ tensor, $I_2$ vector/trator, $I_3$ axial/axitor sectors). Each gives torsion an independent kinetic term.
2. **Curvature-squared invariants** — $R^2$, Ricci², antisymmetric Ricci², Riemann², all evaluated with the *Ricci-Cartan* connection (which mixes torsion in via the Shapiro identity).
3. **Non-minimal couplings to matter** — $R_{[\mu\nu]}F^{\mu\nu}$, parity-odd Chern-Simons-like terms, etc. These are the operators that can *directly* shuttle energy between the gravitational and electromagnetic sectors.

**The total parameter space**: parity-even sector has roughly 10 independent couplings; parity-odd extension ~60. The scan is genuinely high-dimensional.

**The methodological commitment** (worth saying out loud): we don't pick a theory in advance. We let the data — meaning the conversion probability and the causality structure — pick the regions of the Lagrangian that survive.

---

## 7. The question + program

**The question** (read out loud, slowly):
> *Is there any combination of parity-even quadratic curvature and torsion terms whose linearised PDEs, in a background magnetic field, give an amplified Gertsenshtein conversion — and are those surviving theories physically viable?*

**Three things to do for every Lagrangian in the scan:**
1. Derive the linearised field equations symbolically.
2. Screen for causality-violating regimes via all-mode eigenvalue stability.
3. Compute the conversion probability and compare to the GR baseline.

**Pivot line**: "All three steps are what TIDAL automates."

---

## 8. TIDAL — overview

**Acronym, once**: Tensor Integration and Derivation for Any Lagrangian.

**Pipeline** (the tikz diagram on the slide):
TOML → Mathematica (xAct/xTras) → JSON → Python numerics → measurements & inference.

**What the pipeline buys** (in physics-relevant terms):
- A new candidate Lagrangian becomes runnable — *with derived PDEs, causality screening, and a conversion-probability measurement* — within hours, not weeks.
- The symbolic derivation is *exact*; the numerics reproduces the analytic Boccaletti formula to **0.04%**, so deviations from analytic GR can be trusted as physics.

**Core claim** (read verbatim from slide):
> "Every equation the numerics solves is traceable back to a single Lagrangian I typed in. No hand-derived PDEs, no manual index gymnastics."

**Why this matters for *this* project specifically**: scanning ten or sixty couplings systematically is impossible by hand. Without the pipeline, the abstract's "phenomenological scan" reduces to "pick one model and hope." The pipeline turns the abstract's program from a slogan into a procedure.

---

## 9. TIDAL — the symbolic side

**Input**: a TOML file listing fields, Lagrangian, background fields, gauge fixing.

**Processing** (Mathematica + xAct/xTras):
- `xPert` for perturbation around arbitrary backgrounds.
- Euler-Lagrange w.r.t. each field independently (tetrad, spin connection, matter).
- Component decomposition, gauge fixing, mass / coupling matrix extraction.

**Two clever bits worth a sentence each:**

1. **Component-level Euler-Lagrange** (`architecture.tex:372-393`). Decompose to scalar components first, then apply standard calculus E-L per component. Skips the abstract-tensor-index machinery. Without this, the $\widetilde R^2$ derivation (~45 contracted indices) takes 77 minutes; with it, it takes 5 seconds. *That speedup is what makes scanning the parameter space tractable.*

2. **Custom parallelisation framework**. The abstract calls this out specifically. The xAct contraction kernel runs in parallel over independent symbolic sub-expressions. This was necessary because even with component-level E-L, theories with full quadratic torsion + curvature can produce thousands of terms. Without parallelisation, one Lagrangian variant takes a working day; with it, minutes.

**Tagline**: *"Not a lookup table of hand-coded theories — a symbolic compiler for arbitrary Lagrangians."*

---

## 10. The modal solver — architecture and what made it hard

**This is the centrepiece of the numerical side.** Don't dilute it with multi-backend selling — focus on what the modal solver *is* and what we had to fix to make it work for PGT.

**What it does** (one paragraph, say out loud):
For systems linearised around a translation-invariant background, the equations are *linear with constant coefficients in space*. Fourier-transform: each mode $k$ evolves under its own finite-dimensional matrix $M(k)$. Diagonalise $M(k)$ once, compute $y(t,k) = e^{M(k)t}y_0(k)$ for any $t$. **Result: machine-precision accuracy ($|dE/E| \sim 10^{-14}$), no CFL condition, time-evolution cost independent of $t_\text{end}$.**

**Why it's the right tool for *this* problem specifically** (this is the bit to say slowly):
- Every Lagrangian we test linearises to constant-coefficient PDEs in flat space + uniform $B_0$. The modal solver's eligibility criteria are *precisely* the regime of the abstract.
- Conversion probabilities for amplified parameter regions can range over 10+ orders of magnitude (e.g., Stage D1 detected suppression to $A \sim 10^{-12}$). Time-stepping schemes accumulate CFL error that swamps signals below $\sim 10^{-8}$. The modal solver doesn't accumulate any error, ever.
- A nested-sampling inference run calls the solver $\sim 10^5$ times. Modal is $\sim 100\times$ faster than CVODE on this workload (`modal_solver.tex:231-246`). That's the difference between a campaign that finishes and one that doesn't.

**The hard problems we had to solve to make it work for PGT** (this is the war-stories part — pick 2–3 to mention):

1. **Rank-deficient mass matrices.** PGT generically has constraint fields (e.g., $A_0$, redundant torsion gauge DoF). The mass matrix $M$ is singular; you can't invert it. Fix: Schur-complement elimination — split into dynamical and constraint blocks, solve the constraint block algebraically, recover the evolution matrix on the dynamical sector only. Without this, the entire PGT family was unreachable. (`modal_solver.tex:200-222`.)

2. **Doubly rank-deficient blocks** (when constraints couple to constraints). Standard eigendecomposition of $M^{-1}K$ fails. Fix: generalised eigenvalue problem via QZ decomposition (`scipy.linalg.eig(K, M)`).

3. **Hidden block-diagonal structure**. The PGT+EM linearised system splits into independent polarisation blocks (e.g., $h_5\leftrightarrow a_1$ and trace↔torsion sectors). One sector is physical, the other is ghost-unstable. The modal solver detects this and works on the physical sector only.

4. **The Nyquist energy leak** — pseudospectral war story worth telling: eigendecomposition of the highest-frequency Fourier mode produces complex coefficients; the inverse FFT silently discards the imaginary part at the Nyquist bin. Energy was drifting at $1.5\times 10^{-5}$ — small enough to miss casually, six orders of magnitude away from machine precision once you went looking. Zero the Nyquist bin before eigendecomposition: $|dE/E|$ drops to $2\times 10^{-14}$. Standard practice in pseudospectral methods (Boyd 2001) but a rite of passage to re-discover. (`modal_solver.tex:109-144`.)

5. **Non-normal evolution matrices for localised backgrounds** — when $B_0(x)$ is spatially Gaussian rather than uniform, the Fourier-space coupling becomes a convolution, and the resulting $\sim 3000\times 3000$ matrix is *non-normal* (gradient operators $ik$ create eigenvalues with large positive real parts despite conservative physics). Individual $e^{\lambda t}$ overflow even though $e^{At}\cdot y_0$ is bounded. Fix: never form the eigendecomposition; use `expm_multiply` (Al-Mohy & Higham 2011) to compute the action of the matrix exponential directly. 29 s → 1.4 s, and the result agrees with Boccaletti's localised formula to 0.04%.

**Pitch line for this slide**: *"The modal solver gives us an exact answer. Most of the engineering work was making 'exact' actually mean exact in the presence of constraints, gauge symmetries, and the floating-point traps that follow from them."*

---

## 11. Causality and stability — checking the surviving theories

**This is the abstract's "characteristic speeds / causality" slide.** Frame it honestly: the project doesn't (yet) compute characteristic speeds in the formal PDE-theory sense; it does the closely related thing of computing *all-mode eigenvalues of the linearised evolution* and rejecting any region with growing modes.

**The mechanism** (`tidal/inference/_prior_stability.py`, `tidal/measurement/_stability.py`):
1. For every parameter point in the scan, build the modal evolution matrix $M(k)$ at each spatial wavenumber $k$.
2. Compute the largest real part of any eigenvalue, $\max_k\,\mathrm{Re}(\lambda(k))$.
3. If $\max\mathrm{Re}(\lambda) > 0$: the linearised system has a tachyonic / gradient-unstable mode. Reject the parameter point.
4. The remaining region is the *causality-preserving slice* of parameter space.

**Wired into Bayesian inference**: the stability check is part of the prior, not a post-hoc filter. Tachyonic samples get $\log\mathcal{L} = -\infty$, so PolyChord literally never samples there. The cost is $\sim 1$ ms per sample — about 5 s for a 5000-point prior draw. Cheap enough to do at every sample.

**Why this matters for the talk's narrative** (this is the connection to make explicit):
- The very same machinery that detects "this Lagrangian is unstable" also detects "this looked like amplification but was actually tachyonic onset" (issue #238 — see backup).
- Without it, the apparent amplification factor $A\sim 665$ in the early non-minimal sweep would have been claimed as a result. With it, those points are flagged as causality-violating before they enter the posterior.

**The methodological lesson worth landing**: *"In a phenomenological scan, the cheapest sanity check has to be embedded in the prior. Otherwise the inference will spend most of its time in regions where the linearisation is invalid."*

**Future enhancement to mention if asked**: explicit dispersion-relation extraction $\omega(k)$ would give phase/group velocities directly. Listed as Phase I in `NEXT_PHASES.md:405-437`; the eigenvalue stability gate is the current implementation.

---

## 12. Results (placeholder)

- Deliberately blank. Populate after the next sweep round (currently Stage D3 / D2.2 Shapiro / parity-odd survey).
- If forced to present earlier, the safe fallback is the **Boccaletti validation**: numerical agreement to 0.04% with the localised-$B_0$ analytic formula across all $R/\sigma$ ratios. Honest framing — methodology validation, not physics result.

---

## 13. Outlook — *the survey progress map*

The abstract's "systematic phenomenological program" is being executed as a multi-stage HPC campaign. This slide should show the *map* — what's scanned, what's survived, what's queued. **The audience's takeaway should be: "this is a real survey, with results, not a vague aspiration."**

**Stages completed and what they found:**

| Stage | Theory class | Outcome |
|-------|--------------|---------|
| **A** | Dark-photon plasma (T1) | Amplify NULL; **suppress informative** ($D_{\rm KL}=1.98$ nats, $\sim 100\times$ suppression in decoupling corner). |
| **B** | Einstein-Cartan (T2) | NULL — torsion-trace coupling alone is Gertsenshtein-neutral. |
| **C** | $R^2$-PGT structural | Algebraic null in TT channel — no run needed. |
| **D1** | Ricci-EM ($\delta_1 R_{[\mu\nu]}F^{\mu\nu}$, T4) | **Strong suppression: Bayes factor $\sim 10^7$, $A_{\min}\sim 4\times 10^{-12}$.** Mechanism: destructive interference at $|\delta_1|\approx 1.3$. Amplification disfavoured (BF=0.10). |
| **D2.0** | Bahamonde YM-PGT (T5, 5D) | NULL both directions ($\log Z \approx \pm 0.45$ to $\pm 0.62$). YM-PGT sector is Gertsenshtein-neutral. |
| **D2.1** | Barker $\chi$-axial (T5+$\chi$, 6D) | NULL both directions; $\chi$ posterior centred at zero. Adding the Barker coupling changes nothing. |

**The story emerging**: amplification is *elusive*. Multiple natural-looking PGT extensions produce no amplification or *suppression*. The phenomenological scan is *eliminating regions*, which is informative — it constrains where signal-bearing physics could live.

**What's queued / planned:**
- **D2.2 Shapiro PGT** — pending derivation. Extends D2 with additional torsion-curvature cross terms.
- **D2.3 complete-PGT** — pending derivation. Full quadratic Bahamonde-Barker landscape together.
- **D3 parity-odd YM-PGT-CP (T6)** — pending derivation. The parity-odd sector (Chern-Simons-like couplings); the abstract restricts to parity-even but parity-odd is a natural follow-up.
- **Plasma-Proca extensions** — Gertsenshtein in a realistic plasma with effective photon mass; relevant for actual magnetar / IGM observability.
- **T7 Complete-Even-PGT** (~10 couplings) and **T8 Complete-Odd-PGT** (~60 couplings) — contingent on Stage D finding signal in any subsector. Currently no strong amplification signal, so these are gated on D3 / future stages.

**Pitch line**: *"The campaign so far has eliminated the simple PGT extensions as Gertsenshtein amplifiers. The next stages probe the structurally richer corners — parity-odd couplings, complete higher-curvature models, plasma-modified backgrounds. Each is a config-file change in the pipeline, not a rewrite."*

---

## 14. Summary
- Three bullets, then stop. No rambling.
1. Standard Gertsenshtein is astrophysically useless — any usable signal needs new gravitational physics.
2. We're systematically scanning the parity-even quadratic curvature + torsion landscape, with causality screening built into the inference prior.
3. The first half of the scan is done: simple PGT extensions don't amplify; one direction (non-minimal Ricci-EM) gives strong *suppression*; the higher-dimensional, parity-odd, and plasma corners are next.

**Closing line**: *"Whether or not amplification exists in this landscape, we're building the systematic map of where it could live. Thanks — questions and delivery feedback welcome."*

---

## Discovery anecdotes (for sprinkling and Q&A)

### A. The backreaction validity regime — Hwang & Noh's argument

Foundational context for the linearisation. **Hwang & Noh 2023**: backreaction parameter scales as $\sim P/(8\pi^2)$, so $P\ll 1$ is exactly the linearisation-validity condition. The conversion probability *itself* measures when perturbation theory fails.

**Our extension** — the $B_0\to 0$ theorem: we compute the conversion *coefficient* $C_0\equiv\lim_{B_0\to 0} P/B_0^2$. At $B_0=0$ the background is exactly self-consistent (flat Minkowski with zero gauge field is a trivial solution of the *full nonlinear* equations). Backreaction corrections are $O(B_0)$ and vanish in the limit. So $C_0$ is non-perturbative in the background, perturbative in the perturbation amplitude only.

**Pitch line**: *"We're not computing $P$ at some particular $B_0$ and hoping the answer is right — we're computing the coefficient $C_0 = P/B_0^2$ in the limit, where the background is provably exact."*

### B. The naive Hamiltonian is not conserved

A surprising theoretical fact worth a backup-slide or Q&A anecdote.

For systems with constraint fields (gauge theories, PGT, massive gravity), the naive Hamiltonian — the textbook Legendre transform $H = \sum p_i\dot q_i - L$ — is **not conserved**. We measured $\max|dE/E|=1.37$ for massive gravity. Order-unity violation.

**Why it fails**: the Legendre transform assumes all velocities are independent. For constrained systems, constraint-field velocities are algebraically determined by the dynamical fields. The momentum conjugate to a constraint field vanishes by definition ($p_c\equiv \partial L/\partial \dot q_c\approx 0$); the naive $H$ includes terms that are gauge-artefacts, not physical energy.

**The proper treatment**: Dirac-Bergmann constrained Hamiltonian mechanics. Identify primary constraints, derive secondary constraints from consistency, construct the *reduced* Hamiltonian on the physical phase space. Only this reduced $H$ is gauge-invariant and conserved.

**Pitch line**: *"The Legendre transform you learn in your first mechanics course silently breaks for constrained systems. In gauge gravity the constraints are not optional — so this is a routine surprise that shows up in any honest energy diagnostic."*

### C. Tachyonic onset masquerading as amplification (#238)

Best "we thought we had X, turned out to be Y" moment in the project.

- Early sweep of $\delta_1 R_{[\mu\nu]}F^{\mu\nu}$ showed $A\sim 665$ near $\delta_1\approx 1.0$, $\alpha_2\approx -0.91$.
- Time-independence test: $P$ grew $6.5\times$ at $t=20$ and $808\times$ at $t=50$. Exponential, not oscillatory.
- Diagnosis: Schur-complement coupling shift $\Delta m^2_{a_1}\propto\delta_1^2$ drives the photon's effective mass through zero. Tachyonic onset at the boundary of the Hwang validity region.
- Methodological surprise: the $B_0$-scaling sanity check was *blind* to this, because tachyonic growth rates are $B_0$-independent.
- Resolution: the all-mode eigenvalue stability gate (slide 11) catches it before the simulation runs.

### D. Torsion-independence is *algebraic* (#199, #200)

Cleanest physics result we have. The minimal PGT+EM linearised system decouples into two completely independent polarisation channels:
- $h_5\leftrightarrow a_1$ (TT graviton ↔ transverse photon): **torsion-free**, gives standard $\sin^2(\kappa B_0 D/2)$.
- trace sector ↔ $a_2$ ↔ torsion: ghost-unstable, *algebraically inaccessible from TT initial conditions*.

A 200-point scan over all torsion parameters gave $P_\text{max} = 6.249987\times 10^{-6}$ to machine precision — matching the GR baseline.

**Pitch line**: *"Torsion doesn't fail to modify Gertsenshtein for numerical reasons. It fails because the polarisation structure won't let it."*

### E. The $\sqrt{-g}$ subtlety

Subtle GR fact worth one sentence on the symbolic-side slide. When linearising an action you must perturb $\sqrt{-g}\,\mathcal{L}$, not $\mathcal{L}$ alone. The volume element contributes a trace term to the Einstein tensor at second order. Massive gravity goes exponentially unstable if you skip it. We caught the bug because the energy-conservation diagnostic flagged unphysical growth in a theory where the answer was known. (`volume_element_fix.tex`.)

---

## The one non-obvious thing (for the closing summary and Q&A)

This is the sentence I should rehearse most carefully.

> "The phenomenological scan over parity-even quadratic gravity is *eliminating regions*. Most extensions don't amplify Gertsenshtein — many actively suppress it, and several boundaries that look like amplification are tachyonic onset. The map of *where amplification cannot live* is itself a result, because it tells you which theoretical structures the relevant physics has to have."

If I land that in 25 seconds at the end of the talk, the rest of the talk is its setup.
