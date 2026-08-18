# Supervisor Meeting — 8 May 2026

**Period**: 17 April (last meeting) to 8 May 2026

---

## Summary

Key workstreams since the last meeting:

1. **Survey progress and roadmap.** The effective dark-photon model, the minimal Einstein-Cartan theory, the Ricci-EM model, and the entire propagating-torsion nonminimal sector (five nested sub-theories up to the full 9-D joint) are complete. Parity-odd, complete-PGT, and higher-curvature EM remain.
2. **Perturbative reduction — major new research direction.** Handling theories with higher-derivative corrections required developing a novel approach after all standard Hamiltonian methods failed. The equations of motion can now be solved correctly without ghost contamination. For the specific case where the small parameter promotes a non-dynamical field to dynamical, however, the measurement and inference pipeline has not yet been extended — so physics results from that class remain outstanding.
3. **Stability filtering in the linearized regime — open question.** A key methodological decision was made to exclude exponentially-growing modes from the inference. But the physical validity of that choice is an open question.
4. **Practice talk written.** The talk is now ready to schedule and give.

---

## 1. Constraint-torsion theory with a nonminimal Ricci–photon coupling (the only structured result)

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}\tilde{R} + \alpha_i\,I_i + \delta_1\,\tilde R_{[\mu\nu]}\,F^{\mu\nu} - \frac{1}{4}F_{\mu\nu}F^{\mu\nu}$$

The three torsion-mass invariants $\alpha_{1,2,3}$ plus a single nonminimal coupling $\delta_1$ between the Ricci-Cartan tensor and the photon field strength.

**Crucially, this theory has _constraint_ (non-propagating) torsion** — there is no kinetic term for torsion, so it carries no dynamical modes and acts as an auxiliary field. $\delta_1$ is the only term that _connects_ torsion to the Gertsenshtein channel. This is the only completed theory in the _constraint torsion + nonminimal coupling_ quadrant of the survey.

### Results (paired amplify / suppress runs)

| Run      | $\log Z$          | joint $D_{\rm KL}$ | Posterior signal                                                                                                |
| -------- | ----------------- | ------------------ | --------------------------------------------------------------------------------------------------------------- |
| Amplify  | $-2.26 \pm 0.07$  | 1.79 nats          | $A_{\max}=1.26$; Bayes factor 0.10 vs null — **model disfavored 10:1 for amplification**                       |
| Suppress | $+15.92 \pm 0.13$ | 8.91 nats          | $A_{\min} \approx 4\times 10^{-12}$ at MAP; valley reaches $\sim 5\times 10^{-9}$ across $(\alpha_i, \delta_1)$ |

- $\delta_1$'s marginal $D_{\rm KL}$ is **96% of the joint** in amplify — a single coupling carries the entire signal. In suppress, $\delta_1$ leads but $\alpha_{1,2,3}$ each contribute $\sim 0.25$–$0.30$ nats: the deepest suppression valley needs all four parameters to coordinate.
- Suppression depth is far below the analytic estimate ($\sim 10^{-3}$). Consistent with destructive interference; mechanism not yet pinned down.
- **This is the first non-trivial result of the survey** — every other completed theory has been null.

---

## 2. Survey progress and roadmap

### Theory partition

Organizing the completed theories by whether torsion propagates and whether nonminimal couplings are present gives a useful summary of results so far:

|                                      | constraint torsion       | propagating torsion                              |
| ------------------------------------ | ------------------------ | ------------------------------------------------ |
| **minimal (no nonminimal coupling)** | Einstein-Cartan — null   | $R^2$-PGT — $b_5$ decouples structurally         |
| **nonminimal coupling present**      | **Ricci-EM — structure** | YM-PGT family (5 nested sub-theories) — all null |

Plus the phenomenological class (dark-photon-plasma, plasma Gertsenshtein) and not-yet-completed parity-odd / complete-PGT / higher-curvature-EM theories.

**Observation:** the most structure found so far exists in the constraint-torsion case. We cannot yet say whether this reflects something fundamental about the landscape or is simply where we happened to look — propagating-torsion with nonminimal couplings dominates the parameter count of what has been tested, and the parity-odd and complete-PGT sectors remain untested. That said, it is enough to motivate the question of whether the constraint-torsion direction deserves more attention.

Within the YM-PGT family we tested **nested sub-theories** prominent in the literature before the full joint: Bahamonde (arXiv:2402.08937) $\subset$ Barker (arXiv:2406.12826) $\subset$ Shapiro (arXiv:hep-th/0103093) $\subset$ full (9-D). Each layer adds couplings; each layer was null.

### YM-PGT Lagrangian

The full nine-dimensional family is

$$\mathcal{L} = \frac{1}{\kappa^2}\tilde{R} + \beta_i\,I_i - \frac{\xi}{4}(F_T)_{\mu\nu}(F_T)^{\mu\nu} + \delta_1\,\tilde R_{[\mu\nu]}\,F^{\mu\nu} + \chi\,\tilde R_{[\mu\nu]}\,(F_T)^{\mu\nu} + \zeta_a\,(\nabla T)_a^{\mu\nu}\,F_{\mu\nu} - \frac{1}{4}F_{\mu\nu}F^{\mu\nu}$$

where $(F_T)_{\mu\nu} = \partial_\mu T_\nu - \partial_\nu T_\mu$ is the field strength of the torsion-trace vector $T_\mu$, and $(\nabla T)_a^{\mu\nu} F_{\mu\nu}$ collects three independent derivative-torsion–EM contractions.

| Coupling        | Role                                                                                                  |
| --------------- | ----------------------------------------------------------------------------------------------------- |
| $\beta_{1,2,3}$ | _propagating-torsion masses_ (analogue of §1's $\alpha_i$ once $\xi \neq 0$ activates the trace)      |
| $\xi$           | torsion-trace Yang-Mills kinetic term — gives the trace its propagating character                     |
| $\delta_1$      | $\tilde R_{[\mu\nu]} F^{\mu\nu}$ — same nonminimal coupling as in §1                                  |
| $\chi$          | $\tilde R_{[\mu\nu]} (F_T)^{\mu\nu}$ — kinetic mixing of Ricci-Cartan and the torsion field strength  |
| $\zeta_{1,2,3}$ | three independent derivative-torsion–EM cross-terms (Shapiro)                                         |

The nested sub-theories activate progressively: Bahamonde keeps $\beta, \xi, \delta_1$; Barker adds $\chi$; Shapiro adds $\zeta_{1,2,3}$; the full 9-D activates everything.

### Roadmap

| Theory class                                           | Free dim | Status                 | Verdict                                                         |
| ------------------------------------------------------ | -------- | ---------------------- | --------------------------------------------------------------- |
| dark-photon-plasma (effective, Proca torsion)          | 4        | Done                   | Null amplify; genuine suppression at decoupling corner          |
| Einstein-Cartan (minimal constraint torsion)           | 3        | Done                   | Null — torsion structurally decouples                           |
| $R^2$-PGT ($b_5$ minimal-quadratic)                    | 4        | Deferred               | $b_5$ decouples from TT channel structurally                    |
| Ricci-EM (nonminimal, constraint torsion)              | 4        | Done + rerun in flight | **Strong suppressor; structure on $\delta_1$. Rerun deferred.** |
| Bahamonde YM-PGT ($\beta_{1\text{–}3}, \xi, \delta_1$) | 5        | Done                   | Null                                                            |
| Barker YM-PGT ($+\chi$)                                | 6        | Done                   | Null — $\chi$ inert                                             |
| Shapiro YM-PGT ($+\zeta_{1\text{–}3}$)                 | 8        | Done                   | Null — $\zeta_i$ inert                                          |
| full YM-PGT (9-D)                                      | 9        | Done                   | Null — all 6 nonminimal couplings inert                         |
| parity-odd YM-PGT                                      | ~22      | **Next**               | Pending; HPC submission planned                                 |
| complete-even PGT                                      | ~20      | Pending                | Derivation in progress                                          |
| complete-odd PGT                                       | ~30+     | Pending                | Derivation in progress                                          |
| Einstein–Maxwell + higher-curvature EM                 | TBD      | Pending                | Blocked on Wolfram-side xAct issue                              |

### YM-PGT null results (compact)

| Theory     | arXiv          | Amp $\log Z$ | Sup $\log Z$ | Bayes factor | Inert couplings                      |
| ---------- | -------------- | ------------ | ------------ | ------------ | ------------------------------------ |
| Bahamonde  | 2402.08937     | $+0.616$     | $-0.449$     | 2.90         | $\delta_1$                           |
| Barker     | 2406.12826     | $+0.618$     | $-0.447$     | 2.90         | $\delta_1, \chi$                     |
| Shapiro    | hep-th/0103093 | $+0.612$     | $-0.615$     | 3.41         | $\delta_1, \chi, \zeta_{1\text{–}3}$ |
| full (9-D) | —              | $+0.6150$    | $-0.6146$    | 3.42         | $\delta_1, \chi, \zeta_{1\text{–}3}$ |

The propagating-torsion nonminimal sector is **Gertsenshtein-neutral**: no coupling opens a conversion channel. Posterior shape is dominated entirely by the $\beta_{1\text{–}3}, \xi$ stability boundary; all six nonminimal couplings have marginal $D_{\rm KL} < 0.06$ nats.

### Dark-photon-plasma model

The torsion trace vector acts as a Proca dark photon kinetically mixed with the photon, with an effective plasma mass $m_A^2$ on the photon.

| Run      | $\log Z$           | $D_{\rm KL}$ | Verdict                                                                                |
| -------- | ------------------ | ------------ | -------------------------------------------------------------------------------------- |
| Amplify  | $-0.073 \pm 0.007$ | 0.024 nats   | Null amplify                                                                           |
| Suppress | $+0.66 \pm 0.05$   | 1.98 nats    | Genuine suppression at decoupling corner ($m_A^2\approx 0.97, \alpha_3 \approx 0.001$) |

### References

- Barker (arXiv:2406.12826) — Barker-PGT, $\chi$ coupling
- Shapiro (arXiv:hep-th/0103093) — Shapiro derivative couplings $\zeta_{1,2,3}$
- Bahamonde et al. (2024) — propagating-torsion phenomenology programme
- An, Pospelov, Pradler (arXiv:1302.3884) — dark-photon plasma conversion

### Next steps (provisional, subject to discussion)

- **Parity-odd YM-PGT:** tests whether parity-violation in the torsion sector opens what even-parity nonminimal couplings did not.
- **Complete-even and complete-odd PGT:** exhaustive quadratic enumeration. Derivations running locally.
- **Plasma background for the torsion models.** The effective dark-photon model shows that adding an effective photon mass from the plasma environment (Raffelt–Stodolsky resonance) opens conversion channels that are zero in vacuum. The same extension can be applied to the constraint-torsion theories — placing the Ricci-EM (and future constraint-torsion + nonminimal) Lagrangians in a plasma background to ask whether the resonance lifts the suppression or opens amplification channels that the vacuum survey misses.
- **Einstein-Maxwell + Euler-Heisenberg ($F^4$) correction.** The one-loop QED correction to Maxwell adds an $F^4$ term to the Lagrangian. Unlike the $b_5\tilde R^2$ PGT case, this does not promote any constraint field and so the measurement pipeline is unaffected by the Hamiltonian issues. Does the EH nonlinearity modify Gertsenshtein at all, even without any torsion?
- **Possible reorientation toward constraint-torsion + extended nonminimal couplings.** Given the Ricci-EM structure vs the entire propagating-torsion family's emptiness, it may be more fruitful to broaden the _constraint-torsion + nonminimal_ class — single cross-terms in the Ricci-EM mould: scalar-curvature-photon, torsion-trace-to-photon, axial-torsion-to-photon, etc. Worth your input.

---

## 3. Stability filtering in the linearized regime

To run inference on the Ricci-EM class we had to address a fundamental issue: the Padé matrix-exponential that evaluates the field equations is robust for well-conditioned systems, but the eigendecomposition-based path used previously catastrophically rejected all parameter points in the Ricci-EM prior (including physically clean ones) due to ill-conditioning. After switching to a conditioning-robust method, the approach is to compute the growth rate of the solution directly and _exclude_ any parameter point where the linearized fields grow exponentially.

This choice — excluding samples with any exponential growth above a threshold — enabled all the Ricci-EM results above. But it raises a question worth discussing:

**Physical question**: In the linearized regime, an exponentially growing mode will eventually violate the linearization assumption ($\delta g \ll g_{\rm background}$). But should all such modes be discarded? One could argue:

- Large amplification of the Gertsenshtein signal may _require_ some resonant growth to accumulate — the very mechanism by which the photon channel is enhanced might look like tachyonic instability in the linearized equations.
- Some apparent instabilities may be artifacts of the linearized approximation that are cut off in the full nonlinear theory (e.g. by backreaction on the background field, or by nonlinear saturation).
- There may be a middle ground: instabilities that grow slowly enough that the linearized solution remains valid over the physical propagation length of interest, and for which the accumulated conversion is what we actually want to measure.

Currently we treat any growth rate above $\sim 0.15\,\text{s}^{-1}$ as unphysical and return $\log\mathcal{L} = -\infty$. The question is whether this threshold is too conservative, and whether some of the discarded parameter space represents genuine physics.

---

## 4. Perturbative reduction of higher-derivative theories

Some of the most physically interesting PGT Lagrangians include terms that are quadratic in the Riemann-Cartan curvature (e.g. $b_5\tilde R^2$). These generate **fourth-order** equations of motion which, naively, carry Ostrogradsky ghost modes — unphysical negative-energy degrees of freedom arising from the higher-derivative structure.

The standard resolution is _perturbative reduction_: treating $b_5$ as a small coupling, one substitutes the leading-order ($b_5{=}0$) equations of motion back into the correction terms to eliminate the higher time derivatives, leaving a second-order system.

### Methods attempted and why they failed

**(a) JLM substitution (Jaén–Llosa–Molina).**
The most direct approach: substitute the $b_5{=}0$ equations into the $b_5$-correction terms algebraically. This works when the constraint structure of the theory is unchanged at $\mathcal{O}(b_5)$. For our $b_5\tilde R^2$ PGT theory it fails because the correction _promotes constraint fields to dynamical_ — in the $b_5{=}0$ theory, several torsion components are purely algebraically constrained (no time derivatives), but the $b_5$ correction adds kinetic terms for these fields. Substituting the wrong (static) equation of motion for a now-dynamical field gives incorrect correction terms and breaks the perturbative expansion.

**(b) LPS canonical analysis (Lyakhovich–Pluschchay–Sharapov) and Dirac–Bergmann.**
The principled Hamiltonian route: work out the full constraint algebra of the $b_5$ theory and classify all first- and second-class constraints. This is exact but practically infeasible: our 18-field PGT theory with a fourth-order Lagrangian generates an intractable number of constraint equations ($\mathcal{O}(10^3)$ symbolically), and the constraint structure _changes dimensionality_ at $\mathcal{O}(b_5)$ — the phase space itself gains new dimensions when the promoted fields acquire dynamics.

**The key obstruction shared by both:** there is no published recipe in the literature for the _constraint-promotion_ case, where a field that was non-dynamical at leading order becomes dynamical at the next order. Both approaches implicitly assume the number of dynamical degrees of freedom is fixed. This barrier has been recognised as a named open problem in the PGT literature for over 25 years: Blagojević–Nikolić (1983) called it "if-constraints with critical parameters"; Yo–Nester–Ni (1999–2002, gr-qc/9902032 and gr-qc/0112030) called it "constraint bifurcation with field activation"; Blagojević–Cvetković (2018, arXiv:1804.05556) prove explicitly that the perturbative Hamiltonian fails to exist in the relevant subspace. No published method resolves it for the generic PGT case.

### What the equations-of-motion approach gives — and where it stops

For the equations of motion we made progress. Rather than reducing the Lagrangian, we work at the EOM level directly: solve the unperturbed ($b_5{=}0$) equations exactly, then treat the $b_5$-correction as a source $S[y^{(0)}]$ for the second-order base operator $L$. The correction $y^{(1)}$ satisfies $L\,y^{(1)} = S[y^{(0)}]$, solved via the Duhamel convolution integral:

$$y^{(1)}(t) = \int_0^t e^{(t-\tau)A}\,S(\tau)\,y^{(0)}(\tau)\,\mathrm{d}\tau$$

This gives correct, ghost-free equations of motion — the base operator $L$ is always second-order, so no Ostrogradsky modes enter. Verified against exact analytic solutions (Parker–Simon FLRW to $10^{-12}$, driven-oscillator to $10^{-14}$).

**However, solving the equations of motion is not enough to do physics.** Computing observables — conversion probability, energy in each mode — requires a Hamiltonian: you need to know which quadratic combinations of fields and velocities are conserved energy, and how to assign energy to each mode. Without a Hamiltonian, the EOM solution is a set of evolving field amplitudes with no clear physical meaning. For the constraint-promotion case, the Hamiltonian analysis is exactly what fails. The measurement and inference pipeline for these theories is consequently blocked: we can propagate the fields but cannot convert those trajectories to the observables needed for Bayesian inference.

---

## 5. Practice talk

The talk covering this project has been written and is ready to schedule. The structure leads with why the Gertsenshtein effect is astrophysically useless at GR rates, introduces the survey question, and describes the inference framework before the results (yet to be added). Ready to give — looking to schedule a time.

---

## Questions

### 1. Stability filtering — how conservative should we be?

We currently exclude all parameter points where linearized fields grow exponentially above a threshold. This enabled the Ricci-EM inference results. But:

- [ ] Can large amplification of the Gertsenshtein signal occur _without_ some resonant growth? If the conversion mechanism is fundamentally a resonance, the very signal we are looking for might live in the excluded region.
- [ ] Are some of the apparent instabilities artifacts of the linearized approximation, stopped by nonlinear effects in the full theory?
- [ ] Is there a physically motivated threshold below which growth is acceptable — e.g. growth slow enough that amplitudes remain in the linearized regime over the propagation distance?

### 2. Is a Hamiltonian analysis possible at all for the constraint-promotion case?

As described in §4, both JLM substitution and LPS/Dirac–Bergmann fail on the constraint-promotion case ($b_5\tilde R^2$ PGT: h₄, h₇, h₉ becoming dynamical). The failure is not a matter of practical complexity — it is structural. Blagojević and Cvetković (2018) prove it explicitly: "the expressions for $c_n$ when $\bar{b} = 0$ cannot be obtained by taking the limit $\bar{b} \to 0$ of the generic result" — the perturbative Hamiltonian does not exist in this limit.

This has been identified as a named open problem three times in 25 years (Blagojević–Nikolić 1983; Yo–Nester–Ni 1999–2002; Blagojević–Cvetković 2018). We have the Duhamel-source equations of motion (ghost-free, numerically verified to $10^{-12}$), but without a Hamiltonian we have no formal ghost-freedom guarantee beyond leading order and no proper degree-of-freedom count. The most promising lead in the literature is Lyakhovich (2021, arXiv:2102.10579), which has an existence theorem for a generalized Stückelberg recipe — but it has not been applied to PGT curvature-squared theories.

- [ ] Do you know of any framework that handles the case where the number of dynamical degrees of freedom changes order-by-order in perturbation theory? Is there a route — even in principle — to perform a Hamiltonian analysis for this class?

### 3. Campaign direction — constraint torsion vs propagating torsion

The only completed theory with posterior structure has constraint torsion + a nonminimal coupling. The entire propagating-torsion nonminimal sector (five sub-theories) is null.

- [ ] Is it worth pivoting to broaden the constraint-torsion + nonminimal class (Ricci-EM-mould theories with different cross-terms) rather than continuing into propagating-torsion / parity-odd / complete-PGT?
- [ ] Or complete the parity-odd and full-quadratic enumeration first?

### 4. Manuscript appendix — brief introduction to PGT?

PGT isn't covered in the masters curriculum, so the examiner will not have the background context to follow the manuscript. I am wondering whether to include an appendix giving a brief introduction to PGT — collecting key results and a few illustrative figures from canonical references (Blagojević, arXiv:gr-qc/0302040, was particularly helpful for me) into a self-contained overview.

- [ ] Is this the right approach for setting context, or is there a different convention for this?

### 5. HPC queue — is this normal?

Standard-QOS jobs consistently sit PENDING indefinitely (current jobs have been queued for days without running). In practice every job ends up submitted to the INTR queue (1h wall limit), which schedules immediately but times out before convergence, requiring repeated `--resume` resubmissions. The current nested-sampling runs each take 3–5 INTR slots to converge.

- [ ] Is this typical for CSD3 at the moment, or is there something wrong with how the jobs are configured?
- [ ] Is there a recommended way to run jobs that genuinely need more than 1 hour without hitting the standard queue?

### 6. Manuscript §2.2 — TT-gauge-in-matter citation

We evolve the metric perturbation without an equation-of-motion gauge constraint, and impose only TT-compatible initial conditions in the vacuum region where the wavepacket originates. Per issue #167 this is the supervisor-guided approach, but the manuscript prose currently has no published reference for it.

- [ ] Is there a published reference for this TT-in-matter approach (notes? a recent paper?), or should we cite it as "private communication"?

### 7. Carry-forward from 17 April

- [ ] **Scheduling the practice talk** — when?
- [ ] **Sven Krippendorf** — any reply since you reached out?
- [ ] **Will Handley** — anything new on his upcoming PhD positions, and the right timing for an application?
