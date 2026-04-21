# Speaker notes — practice talk (2026-04)

Raw research material for the spoken script. One heading per slide. Phase D turns these into `\note{}` blocks in the Beamer source; Phase E turns them into continuous prose.

All citations are `file:line` references into the codebase so I can dig in later.

---

## 1. Title
- Introduce yourself + project title in one sentence.
- Name the practice framing: "not the viva, rehearsal for feedback."

## 2. Motivation — why this room should care
- GWs are observable; next frontier = high-frequency + BSM couplings.
- The one non-obvious thing I want to land (the payoff for the whole talk):
  > "You cannot simply *add torsion* to the Gertsenshtein effect — the theory itself resists it. That resistance has a structural reason I can show you."
- Three-beat promise: (1) the channel, (2) the question, (3) the tool I built to answer it.

---

## 3. Gertsenshtein — the effect

**Setup to say out loud:**
- Einstein-Maxwell linearised around flat Minkowski; uniform background $B_0$ along $x$ (via Coulomb gauge $\bar A_y = -B_0 z$).
- Expand metric perturbation $h$ (graviton) and EM perturbation $a$ (photon) to first order.
- TT gauge on $h$, Lorenz on $a$, plane-wave reduction $\partial_x = \partial_y = 0$.

**The coupled system** (from `docs/tex/gertsenshtein.tex`):
$$\partial_t^2 h_+ = \partial_z^2 h_+ - \tfrac12 \kappa^2 B_0^2 h_+ - \kappa^2 B_0\,\partial_z a_y$$
$$\partial_t^2 a_y = \partial_z^2 a_y + B_0\,\partial_z h_+$$
- Asymmetric coupling: gravity kinetic normalisation $1/(2\kappa^2)$ vs Maxwell $1$.
- Beat frequency $\Delta\omega = \kappa B_0$ is the *geometric mean* of the two asymmetric couplings — worth mentioning as a sanity check.

**Headline formula**: $P(g\to\gamma) = \sin^2(\tfrac12 \kappa B_0 D)$ with $D = ct$.

**The physics intuition** (the analogy to use): same math as two-state QM mixing, same as axion-photon Primakoff. Background $B_0$ is a "mixing medium" that detunes the graviton and photon dispersion just enough to drive coherent oscillation.

**Deeper, if I have time**: the underlying reason is that the background $B_0$ contributes to the *graviton's* effective dispersion (through its stress-energy), while the photon stays massless. The detuning sets the beat frequency. This framing is more physical than "two-state mixing" and worth one sentence.

---

## 4. Literature + why it's topical now

**Timeline to walk through:**
- **Gertsenshtein 1962** — the original paper, identified the effect and gave the conversion formula.
- **Boccaletti, De Sabbata, Fortini, Gualdi 1970** — rederivation with plasma corrections; this is the formula most modern papers actually use.
- **Raffelt & Stodolsky 1988** — unified with axion-photon mixing in a single mixing-matrix framework. This is how particle physicists cite it.
- **2020s revival** — Aggarwal et al. *Living Reviews in Relativity* 2021 HFGW review; Berlin, Blas, D'Agnolo, Ellis, Harnik, Kahn, Schutz 2022 (microwave cavity detectors); Domcke & Garcia-Cely 2021 (radio telescopes as HFGW detectors).

**Why it's topical in 2026** (three bullets on the slide, say out loud):
1. Astrophysical $B$-fields — magnetars $\sim 10^{15}$ G, intergalactic medium, early universe — are natural labs.
2. It's the *only* linear-in-$h$ EM coupling, so any HFGW cavity-detector scheme is using Gertsenshtein (or its inverse).
3. Model-independent at GR level, so it's a clean place to test *model-dependent* BSM signatures.

**The literature-error anecdote worth dropping if I have 20 seconds** (`docs/tex/gertsenshtein_formula.tex:26-95`):
> Palessandro & Rothman 2023 quote the formula with a different prefactor — off by $\sqrt{4\pi}$ — because they use a non-canonical graviton kinetic normalisation. Dandoy, Lella et al. 2024 (arXiv:2406.17853) confirm our formula independently. The fact that this effect is *still* being rederived in 2024 and people *still* disagree on factors of $\sqrt{4\pi}$ tells you it's under-explored.

---

## 5-6. Torsion via Poincaré gauging

**Gauge idea, for the audience:**
- GR uses the Levi-Civita connection (metric compatible, torsion-free). Why? Historical choice, not a consequence.
- Utiyama 1956, Kibble 1961, Sciama 1962: gravity as a gauge theory of the Poincaré group, analogous to gauging $U(1)$ for QED.
- Gauging the two subgroups:
  - Lorentz → spin connection $\omega^{ab}{}_\mu$ → curvature $R^{ab}{}_{\mu\nu}$.
  - Translations → tetrad $e^a{}_\mu$ → **torsion** $T^a{}_{\mu\nu}$.
- Tagline: *"GR keeps one of the two field strengths of the Poincaré group. PGT keeps both."*

**Geometric intuition** (slide 6):
- Parallel-transport around a closed loop has two independent mismatches: rotational (curvature) and translational (loop fails to close — torsion).

**Why take PGT seriously** (four quick bullets):
1. Couples naturally to fermion spin (Einstein-Cartan, Hehl et al. 1976).
2. Minimal extension — smallest well-motivated deviation from GR.
3. High-curvature / early-universe regimes where GR might fail first.
4. Decades of literature: Sezgin & van Nieuwenhuizen 1980, Nikiforova 2009, Barker 2024 — all identifying ghost-free sectors.

**Honest caveat bullet**: PGT Lagrangians generically have ghost / tachyon sectors. This isn't a footnote; it will drive the results story later. Name-drop Ostrogradsky's theorem if I'm feeling confident.

---

## 7. The question

**Slow down here, read the blockquote out loud:**
> *In a PGT + Einstein-Maxwell Lagrangian, does torsion modify the Gertsenshtein conversion probability, and if so, through which coupling and in which regime?*

**Sub-question** (the more interesting framing): *what does the Lagrangian have to look like for torsion to matter at all?*

**Pivot line**: "To answer this at the level of an arbitrary Lagrangian — and to do it often enough to sweep parameter space — I built TIDAL."

---

## 8. TIDAL — overview (start of the centrepiece)

**Expand the acronym once, never again**: Tensor Integration and Derivation for Any Lagrangian.

**Pipeline** (the tikz already on the slide): TOML → Mathematica (xAct/xTras) → JSON → Python numerics → sweeps/inference/plots.

**What the pipeline buys you — pitch in terms of physics capability, not code metrics:**
- I can write down a new Lagrangian in a config file and have linearised PDEs, a numerical simulation, and a parameter sweep running within hours — not weeks. For the Gertsenshtein investigation, this meant testing dozens of torsion coupling variants systematically instead of committing to one hand-derived model and hoping it was the right one.
- The symbolic derivation is *exact* — no truncation, no hand-simplification errors. The numerical solver can reproduce known analytic results (e.g. the Boccaletti Gertsenshtein formula) to **0.04% precision**, so when I see a deviation from the analytic answer, I can trust that it's physics, not numerics.

**The core claim** to repeat verbatim from the slide:
> "Every equation the numerics solves is traceable back to a single Lagrangian I typed in. No hand-derived PDEs. No manual index gymnastics."

**Why it matters for this project**: varying the Lagrangian is the whole scientific method here. I am going to vary it dozens of times — different torsion operators, different couplings, different gauges, plasma backgrounds. Hand-deriving each variation is intractable and error-prone. And as we'll see, several of the key findings in this project came from being *able* to vary the theory systematically rather than committing to a single model upfront.

---

## 9. TIDAL — symbolic side

**The four cleverest bits** (from `docs/tex/architecture.tex` and `docs/tex/pipeline.tex`):

1. **Component-level Euler-Lagrange** (`architecture.tex:372-393`). The big win:
   > "Decompose the Lagrangian into scalar components *first*, then apply standard calculus E-L per component. Skip the abstract-tensor-index machinery entirely."
   > **900× speedup** over abstract-index VarD for higher-curvature theories (e.g., the $\widetilde R^2$ term with ~45 contracted indices: **77 minutes → 5 seconds** for 2+1D).

2. **Deferred field canonicalisation** (`pipeline.tex:46-146`). Multi-field theories with derived tensors (Maxwell $F_{\mu\nu}$, torsion field strength) stay abstract all the way through `xPert`, expand only *after* perturbation is complete. Fixed three subtle pathologies in one change: #218 photon EOM merge, #250 graviton $h_5$ fragment, #255 torsion cross-sector interference. Collapsed three competing workarounds (~175 LoC) into one unified per-component canonicalisation.

3. **Volume element fix** (`volume_element_fix.tex`). Linearise the action via *native* $\sqrt{-g}$ perturbation — not just $\delta^2 \mathcal{L}$ but $\delta^2(\sqrt{-g}\,\mathcal{L})$. Massive gravity was **exponentially unstable** before this, because the Einstein tensor was missing its trace contribution. Fix: use xPert's native `Perturbation[Sqrt[-Detg[]], n]` support (commit `ca256cb`). **This is a genuinely good anecdote — a subtle symbolic bug that manifested as unphysical growth.**

4. **Ostrogradsky reduction** (`architecture.tex:395-442`). Higher-derivative Lagrangians get auto-reduced to 2nd-order systems via auxiliary fields. Example: graviton-torsion (37 fields after 3 Ostrogradsky auxiliaries) with a singular mass matrix handled by a three-level elimination (mass eigendecomposition → jerk substitution → Schur constraint elimination).

**Processing summary** (one sentence for the slide):
- `xPert` perturbation around arbitrary backgrounds.
- Euler-Lagrange w.r.t. each field independently (tetrad, spin connection, matter fields).
- Component decomposition, gauge fixing, mass / coupling matrix extraction.

**Tagline**: *"Not a lookup table of hand-coded theories — a symbolic compiler."*

---

## 10. TIDAL — numerical side

**The modal solver is the hero.** (`docs/tex/modal_solver.tex`). Four facts worth citing out loud:

1. **Machine precision** ($\sim 10^{-14}$ error from eigendecomposition alone).
2. **No CFL condition.** Exact solution for any $t$ via $\exp(A\cdot t)\cdot y_0$ — cost is $O(1)$ in simulation time. Long-time runs don't cost more than short ones.
3. **Performance numbers worth saying aloud:**
   - Coupled scalars ($N=256$): **1,451× faster** than CVODE (0.003s vs 4.27s).
   - Gertsenshtein ($N=512$, gradient coupling): **83-86× faster**.
   - Localised Gertsenshtein (position-dependent coefficients): **86× faster** via sparse convolution matrix.
   - Proca 3D with constraints: **238× faster**.
   - Long runs ($t_\text{end}=500$): **7.7× faster** — modal is constant, CVODE scales linearly with $t$.
4. **Cross-check capability**: "I can run the same physics through IDA, CVODE, leapfrog, and the modal solver and have them agree to machine precision. When they disagree, I know there's a bug — and I've caught several that way."

**The Nyquist mode bug** (`modal_solver.tex:109-144`) — one of the best bug stories in the project:
> Eigendecomposition of the highest-frequency Fourier mode creates complex coefficients. `irfft` silently discards the imaginary part at the Nyquist bin (it must be real for a real signal). Energy was drifting at $1.5 \times 10^{-5}$ — small enough to miss casually, but a clear fingerprint once you went looking. **Zero the Nyquist bin before eigendecomposition and you're at $2\times 10^{-14}$ — machine precision.** This is standard practice in pseudospectral methods (Boyd 2001) but is a rite of passage to re-discover on your own.
> **Worth telling verbatim as a backup-slide anecdote.**

**Rank-deficient mass matrix / spurious tachyons** (`modal_solver.tex:582-697`, issue #256):
> PGT with $\widetilde R^2$ produces a singular mass matrix (hidden constraints from Poincaré symmetry). The old modal code mis-classified acceleration operators and spawned fake tachyonic modes with $\mathrm{Re}(\lambda) \approx 1$ that exploded. Fixed in #256 via a unified three-level elimination using QZ decomposition.
> This one is thematically important for the talk: *the investigation was literally blocked by a modal-solver artefact that masqueraded as physics*. Worth a backup slide.

**Analytical Jacobian tiers** (if there's time): dense ($N\leq 2k$, 5.3×), sparse ($2k<N\leq 200k$, 2.5×), GMRES ($N>200k$).

---

## 11. TIDAL — sweeps, inference, HPC

**Declarative downstream pipeline:**
- `tidal sweep` — 13 measurement types, grid/LHS/Sobol/Morris.
- `tidal sample` — Bayesian inference via MC or nested sampling (dynesty dev / PolyChord HPC), with priors (uniform, log-uniform, normal, arctan-uniform), hard constraints, and three likelihoods (maximise, gaussian, threshold).

**HPC numbers worth citing** (`CLAUDE.md`, HPC section):
- CSD3 (Cambridge), Sapphire (112 cores) / Icelake / CClake partitions.
- Sweep parallelism sweet spot: **~98% of ideal efficiency** at `--parallel 32` on a 90-point plasma-Gertsenshtein benchmark.
- Super-linear speedup at $P\in\{8,16\}$ from BLAS cache locality.
- Real campaign number: the null-result sweep was **276 runs** across $(\xi, \alpha, \delta_m)$ with $10^{-6}$ precision agreement with the analytic formula.

**Land the block-quote line** verbatim:
> "The same Lagrangian that gave me a single simulation can give me a 276-point, five-parameter sweep with Bayesian uncertainty quantification — without writing new code."

---

## 12. Results (placeholder)

- Deliberately blank. Populate after the next sweep round: plasma-Proca, dark photon, full PGT HPC.
- If forced to present earlier, the safe fallback is the **Boccaletti validation**: our numerical solver agrees with the localised-$B_0$ Gertsenshtein formula to **0.04% mean error** across all $R/\sigma$ ratios (`gertsenshtein.tex:161-200`) — and to **<0.015 RMS** across a 40-point $B_0$ sweep at $N=1024$.
- This is "numerical methodology validation" framing, not "physics result" framing. Honest.

---

## 13. Outlook — what's coming

Three concrete upcoming investigations:
1. **Plasma-Proca Gertsenshtein** — realistic plasma background (effective photon mass); magnetars, neutron-star magnetospheres, IGM.
2. **Dark-photon BSM coupling** — kinetic mixing with a hidden $U(1)$; cross-references HFGW dark-photon searches.
3. **Full PGT Lagrangian HPC sweep** — parameter-space map of torsion couplings on CSD3.

Payoff line: *"Each of these is now a config-file change, not a rewrite."*

---

## 14. Summary
- Three bullets, then stop. No rambling.
- Gertsenshtein is old, under-explored, topical again.
- Torsion is the half of Poincaré-gauge gravity GR discards; asking whether it modifies Gertsenshtein is the natural first BSM question.
- TIDAL makes the question tractable, and the physics answers are coming.

---

## Physics discoveries and insights worth telling

Raw material for anecdotes to sprinkle into slides or pull into Q&A. Focus is on what's *physically surprising* — things an astrophysicist would find interesting about the interplay between theory, linearisation validity, and constrained systems.

### 1. When does linearisation break down? The backreaction regime and the Hwang argument

This is foundational context for the whole project and worth a sentence on the Gertsenshtein intro slides or the TIDAL slides.

**The problem**: we linearise around flat Minkowski + uniform $B_0$. But a uniform $B_0$ carries stress-energy that should curve spacetime — the background isn't self-consistent. When is the linearisation trustworthy?

**Hwang & Noh (2023)** gave the clean answer: the backreaction parameter scales as $\sim P/(8\pi^2)$. Backreaction is negligible whenever the conversion probability $P \ll 1$ — the conversion probability *itself* measures when perturbation theory fails. For realistic magnetar fields ($P \sim 10^{-10}$), backreaction is utterly negligible.

**Our extension — the $B_0 \to 0$ theorem**: The conversion coefficient $C_0 \equiv \lim_{B_0 \to 0} P/B_0^2$ is computed on an *exactly self-consistent* background. At $B_0 = 0$, all field strengths vanish, so flat Minkowski with zero gauge field and zero torsion is a trivial solution of the *full nonlinear* equations. All backreaction and EFT corrections enter at $O(B_0)$ or higher and vanish in the limit. Since graviton-photon coupling requires $B_0$ (angular momentum conservation), $P \propto B_0^2$ analytically, giving a finite $C_0$ that is exact.

**Pitch-line**: *"We're not computing $P$ at some particular $B_0$ and hoping the answer is right. We're computing the coefficient $C_0 = P/B_0^2$ in the limit $B_0 \to 0$, where the background is provably exact. The result is non-perturbative in everything except the perturbation amplitude."*

### 2. Tachyonic onset masquerading as amplification (#238) — and the lesson about validity

The most instructive physics discovery in the project. The story:

- Swept the non-minimal coupling $\delta_1 \widetilde R_{[\mu\nu]}F^{\mu\nu}$; measured amplification factor $A = P_\text{torsion}/P_\text{GR} \sim 665$ near $\delta_1 = 1.0$, $\alpha_2 \approx -0.91$.
- This looked like a strong signal. But a time-independence test revealed exponential growth: $P$ grew from $6.5\times$ at $t=20$ to $808\times$ at $t=50$. **Exponential, not oscillatory.**
- **Physics diagnosis**: the non-minimal coupling introduces a Schur-complement shift to the photon's effective mass: $\Delta m^2_{a_1} \propto \delta_1^2$. At certain parameter values, $m^2_\text{eff}(a_1)$ passes through zero. When $m^2_\text{eff} < 0$, the photon mode is **tachyonic** — growth $\sim e^{|\gamma|t}$, not oscillation.
- **Why it connects to Hwang**: this is the tachyonic instability Hwang & Noh identified as the fundamental limit of the linearisation. When $P \gtrsim 1$, the perturbation is large enough to backreact and the linearised treatment breaks down. But torsion couplings can drive $P \to 1$ at *much smaller* $B_0$ than GR does, because the tachyonic growth rate is $B_0$-independent — it's set by the effective mass, not the coupling strength.
- **Methodological surprise**: the $B_0$-scaling test (measuring $C_0 = P/B_0^2 = \mathrm{const}$) *cannot* distinguish amplification from instability, because the tachyonic growth rate doesn't depend on $B_0$. A sanity check we trusted turned out to be blind to exactly the artefact it was supposed to catch.
- **Pitch-line**: *"We had to go beyond the standard toolkit of checks to detect this — it required a time-independence test at multiple $t$ values, which is not standard practice in this literature."*
- Source: `docs/AMPLIFICATION_INVESTIGATION.md:162-196`.

### 3. The naive Hamiltonian is not conserved — and why that's physically surprising

This one is a good backup-slide or Q&A anecdote that would genuinely surprise the audience.

**The claim**: for systems with constraint fields (gauge theories, PGT with auxiliary torsion modes, massive gravity), the naive Hamiltonian — the one you get by performing a straightforward Legendre transform $H = \sum_i p_i \dot q_i - L$ — is **not a conserved quantity**. We measured $\max|dE/E| = 1.37$ for massive gravity. Order-unity violation.

**Why it fails**: the Legendre transform assumes all velocities are independent. For constrained systems, constraint-field velocities are algebraically determined by the dynamical fields, not independent degrees of freedom. The momentum conjugate to a constraint field vanishes by definition ($p_c \equiv \partial L / \partial \dot q_c \approx 0$), so the naive $H$ includes terms that are gauge-artefacts, not physical energy.

**The proper treatment**: Dirac-Bergmann constrained Hamiltonian mechanics. You must identify primary constraints ($p_c \approx 0$), derive secondary constraints from consistency ($\dot p_c \approx 0$ gives Gauss's law for EM, for example), and construct the *reduced* Hamiltonian on the physical phase space. Only this reduced $H$ is gauge-invariant and conserved.

**Why astrophysicists should care**: any gravitational theory with extra fields (PGT, Horndeski, massive gravity, etc.) will generically have constraint sectors. If you naively compute $H$ from the full Lagrangian and check energy conservation, you may conclude the simulation is broken when the physics is actually correct — or worse, conclude everything is fine when it isn't.

**Our workaround**: restrict the Hamiltonian to the dynamical sector only (graviton + photon terms), excluding constraint fields. This gives the correct conversion probabilities because the observable is *relative* energy transfer between sectors, and constraint contributions cancel.

**Pitch-line**: *"The Legendre transform that you learn in your first mechanics course — the one that gives you $H$ from $L$ — silently breaks for constrained systems. And in Poincaré Gauge Theory, the constraints are not optional; they're built into the gauge structure."*
- Source: `constraint_fields.tex:373-475`, issue #178.

### 4. Torsion-independence is *algebraic*, not numerical (#199, #200)

Worth embedding in the results or backup slides. This is the cleanest physics result we have.

- The minimal PGT + EM system decomposes into two completely decoupled channels after plane-wave reduction:
  - $h_5 \leftrightarrow a_1$: torsion-free, stable, gives the standard $\sin^2(\kappa B_0 D/2)$.
  - trace sector $\leftrightarrow a_2 \leftrightarrow$ torsion: torsion-dependent, ghost-unstable, algebraically inaccessible from TT initial conditions.
- This is *exact* and *algebraic*, not numerical. A 200-point sweep over all torsion parameters gave $P_\text{max} = 6.249987\times 10^{-6}$ to machine precision, matching the GR baseline.
- **The physical reason**: the TT graviton ($h_+$) and the transverse photon ($a_x$) live in a polarisation sector that has *zero overlap* with any torsion mode. Torsion couples only to the trace/longitudinal sector, which a vacuum GW source doesn't excite.
- **Pitch-line**: *"Torsion doesn't fail to modify Gertsenshtein for numerical reasons. It fails because the polarisation structure won't let it — TT gravitons simply don't talk to torsion modes."*

### 5. The propagating-torsion dead end: coupling and instability from the same operator (#236)

Quadratic PGT with propagating torsion ($\xi F_T^2$ kinetic term) plus $\delta_1 \widetilde R_{[\mu\nu]}F^{\mu\nu}$:
- 5D parameter scan: 497/500 points unstable. Stability window: $|\delta_1| < 0.005$.
- Within the stability window: $A = 1.000$ exactly — zero amplification.
- **The physics**: $\widetilde R^2$ is what makes torsion propagate (gives it a kinetic term), but it's *also* the source of Ostrogradsky ghosts (4th-order time derivatives). The coupling that would let torsion modify Gertsenshtein is bundled with the instability that destroys the regime you'd measure it in. This is a consequence of Ostrogradsky's theorem, not a numerical accident.
- **Literature context**: Sezgin & van Nieuwenhuizen (1980), Nikiforova (2009), Barker (2024) all identify *sector-specific* ghost-free windows, never a universal closed-form condition. The theory space is genuinely constrained.
- **Pitch-line**: *"Within quadratic PGT, the operator that would do something is also the operator that destroys stability. You can't have one without the other."*

### 6. Volume element linearisation — a subtle point about $\sqrt{-g}$

Worth one sentence on the symbolic-side slide or as a backup-slide detail.

When linearising an action, you must perturb $\sqrt{-g}\,\mathcal{L}$, not just $\mathcal{L}$. The volume element $\sqrt{-g}$ itself depends on the metric perturbation and contributes a trace term to the Einstein tensor at second order. If you omit it (linearise $\mathcal{L}$ alone), massive gravity becomes exponentially unstable. This is a standard GR subtlety but easy to miss in a symbolic pipeline — we caught it because the energy conservation diagnostic flagged unphysical growth in an otherwise well-understood theory.
- Source: `volume_element_fix.tex`.

---

## The one non-obvious thing (for the final summary and Q&A)

**Pulled from the two explore reports, this is the deepest thing I learned from all this documentation:**

The Gertsenshtein effect is a *generic* consequence of the graviton acquiring an effective mass-like dispersion from the background EM stress-energy, while the photon stays massless. The detuning drives coherent oscillation.

Torsion theories *add* new fields and new coupling paths — but **all minimal PGT constructions are blocked from modifying this channel by algebraic structure** (the polarisation block-diagonal decomposition). To break that structure you need non-minimal operators like $\widetilde R_{[\mu\nu]} F^{\mu\nu}$. But those operators generically introduce tachyonic instabilities at their own coupling boundaries.

There is a *structural tension*: genuine propagating torsion (via higher-derivative terms) *always* comes with instabilities in PGT, because ghost-freedom is so constraining. Sezgin & van Nieuwenhuizen 1980, Nikiforova 2009, Barker 2024 all identify sector-specific ghost-free windows, never a universal closed-form condition.

**The clarifying takeaway for astrophysicists in the audience:**
> "You cannot simply 'add torsion' to the Gertsenshtein effect. The theory resists it. Either the polarisation structure blocks the coupling, or the coupling comes bundled with an instability that destroys the regime you wanted to probe. Mapping where that resistance gives way — and what's left of the physics on the other side — is what this project is about."

This is the sentence I should rehearse most carefully. If I can land it in 25 seconds at the end of the talk, the rest of the talk is a setup for it.
