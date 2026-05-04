# Supervisor Meeting — 8 May 2026

**Period**: 17 April (last meeting) to 8 May 2026

---

## Summary

Five workstreams since the last meeting:

1. **T4 Ricci-EM is the live wire.** First completed theory with posterior structure. Initial runs strongly disfavour amplification (Bayes factor 0.10); a rerun with corrected initial conditions appears to amplify ($\Delta \log Z = +4.4$ nats). Verdict held until the suppress cross-check lands.
2. **Survey progress and roadmap.** The effective dark-photon model, the minimal Einstein-Cartan theory, the Ricci-EM model, and the entire propagating-torsion nonminimal sector (five nested sub-theories up to the full 9-D joint) are complete. Parity-odd, complete-PGT, and higher-curvature EM remain.
3. **Stability filtering in the linearised regime — open question.** A key methodological decision was made to exclude exponentially-growing modes from the inference; this enabled results on the Ricci-EM class. But the physical validity of that choice is an open question worth discussing.
4. **Perturbative reduction — major new research direction.** Handling theories with higher-derivative corrections required developing a novel approach after all standard Hamiltonian methods failed. This is now working numerically, with one unresolved theoretical gap.
5. **Practice talk written.** The talk is now ready to schedule and give.

---

## 1. T4 Ricci-EM (the live wire)

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}R + \alpha_i\,I_i + \delta_1\,\tilde R_{\mu\nu}\,F^{\mu\nu} - \frac{1}{4}F_{\mu\nu}F^{\mu\nu}$$

The three torsion-mass invariants $\alpha_{1,2,3}$ plus a single nonminimal coupling $\delta_1$ between the Ricci-Cartan tensor and the photon field strength.

**Crucially, this theory has *constraint* (non-propagating) torsion** — there is no kinetic term for torsion, so it carries no dynamical modes and acts as an auxiliary field. $\delta_1$ is the only term that *connects* torsion to the Gertsenshtein channel; without it, torsion is structurally decoupled (the takeaway from the 17 April meeting). This makes T4 the only completed theory in the *constraint torsion + nonminimal coupling* quadrant of the survey — see §2.

### Results (paired amplify / suppress runs)

| Run | $\log Z$ | joint $D_{\rm KL}$ | Posterior signal |
|---|---|---|---|
| Amplify | $-2.26 \pm 0.07$ | 1.79 nats | $A_{\max}=1.26$; Bayes factor 0.10 vs null — **model disfavoured 10:1 for amplification** |
| Suppress | $+15.92 \pm 0.13$ | 8.91 nats | $A_{\min} \approx 4\times 10^{-12}$ at MAP; valley reaches $\sim 5\times 10^{-9}$ across $(\alpha_i, \delta_1)$ |

- $\delta_1$'s marginal $D_{\rm KL}$ is **96% of the joint** in amplify — a single coupling carries the entire signal. In suppress, $\delta_1$ leads but $\alpha_{1,2,3}$ each contribute $\sim 0.25$–$0.30$ nats: the deepest suppression valley needs all four parameters to coordinate.
- Suppression depth is far below the analytic estimate ($\sim 10^{-3}$). Consistent with destructive interference; mechanism not yet pinned down.
- **This is the first non-trivial result of the survey** — every other completed theory has been null.

### Corrected-IC rerun (in flight, deferred)

Same Lagrangian, but with the initial condition wavevector aligned to the canonical Fourier-grid mode ($k_{\rm IC} = 2\pi/L$, $L=100$) and a different stability criterion.

- Amplify rerun: $\log Z = +2.135 \pm 0.059$, a **$+4.4$ nat shift**. New MAP at $\delta_1 \approx +1.94$, apparent $A_{\max} \approx 42$.
- Caveat: $\sim$26% of the posterior has $P_{\max} > 0.3$, so perturbative validity is borderline at the MAP.
- Suppress cross-check still running. **Verdict held until it lands.**

If the conclusion really does flip with initial-condition choice, that itself is a methodological finding and a question for the meeting (see §Questions).

### Literature

| Claim | Citation | Reference |
|---|---|---|
| Ricci $\times F$ as the sole single-coupling channel | This work — informed by 17 April meeting | T4 derivation |
| Constraint vs propagating torsion in PGT | Hehl et al. 1995 | _Phys. Rep._ 258, 1 |

---

## 2. Survey progress and roadmap

### Theory partition

The linearised landscape splits naturally on two axes — *whether torsion propagates* and *whether nonminimal couplings are present*:

| | constraint torsion | propagating torsion |
|---|---|---|
| **minimal (no nonminimal coupling)** | Einstein-Cartan — null | $R^2$-PGT — $b_5$ decouples structurally |
| **nonminimal coupling present** | **T4 Ricci-EM — structure** | YM-PGT family (5 nested sub-theories) — all null |

Plus the phenomenological class (dark-photon-plasma, plasma Gertsenshtein) and not-yet-completed parity-odd / complete-PGT / higher-curvature-EM theories.

**Emerging pattern:** the only quadrant with structure is *constraint torsion + nonminimal coupling*. The entire propagating-torsion nonminimal sector (five sub-theories, six couplings: $\delta_1, \chi, \zeta_{1,2,3}$, plus kinetic $\xi$) is inert. The kinetic structure that gives torsion its own modes appears to dilute or close the cross-channel — suggesting the constraint-torsion direction may be more fruitful.

Within the YM-PGT family we tested **nested sub-theories** prominent in the literature before the full joint: Bahamonde $\subset$ Barker $\subset$ Shapiro $\subset$ full. Each layer adds couplings; each layer was null.

### Roadmap

| Theory class | Free dim | Status | Verdict |
|---|---|---|---|
| T1 dark-photon-plasma (effective) | 4 | Done | Null amplify; genuine suppression at decoupling corner |
| T2 Einstein-Cartan (minimal constraint torsion) | 3 | Done | Null — torsion structurally decouples |
| T3 $R^2$-PGT ($b_5$ minimal-quadratic) | 4 | Deferred | $b_5$ decouples from TT channel structurally |
| T4 Ricci-EM (nonminimal, constraint torsion) | 4 | Done + rerun in flight | **Strong suppressor; structure on $\delta_1$. Rerun deferred.** |
| T5-Bahamonde ($\beta_{1\text{–}3}, \xi, \delta_1$) | 5 | Done | Null |
| T5-Barker ($+\chi$) | 6 | Done | Null — $\chi$ inert |
| T5-Shapiro ($+\zeta_{1\text{–}3}$) | 8 | Done | Null — $\zeta_i$ inert |
| T5 full (9-D) | 9 | Done | Null — all 6 nonminimal couplings inert |
| T6 YM-PGT-CP (parity-odd) | $\sim$22 | **Next** | Pending; HPC submission planned |
| T7 Complete-Even-PGT | $\sim$20 | Pending | Derivation in progress |
| T8 Complete-Odd-PGT | $\sim$30+ | Pending | Derivation in progress |
| Einstein–Maxwell + higher-curvature EM | TBD | Pending | Blocked on Wolfram-side xAct issue |

### YM-PGT null results (compact)

| Theory | Amp $\log Z$ | Sup $\log Z$ | Bayes factor | Inert couplings |
|---|---|---|---|---|
| T5-Bahamonde | $+0.616$ | $-0.449$ | 2.90 | $\delta_1$ |
| T5-Barker | $+0.618$ | $-0.447$ | 2.90 | $\delta_1, \chi$ |
| T5-Shapiro | $+0.612$ | $-0.615$ | 3.41 | $\delta_1, \chi, \zeta_{1\text{–}3}$ |
| T5 full | $+0.6150$ | $-0.6146$ | 3.42 | $\delta_1, \chi, \zeta_{1\text{–}3}$ |

The propagating-torsion nonminimal sector is **Gertsenshtein-neutral**: no coupling opens a conversion channel. Posterior shape is dominated entirely by the $\beta_{1\text{–}3}, \xi$ stability boundary; all six nonminimal couplings have marginal $D_{\rm KL} < 0.06$ nats.

### Dark-photon-plasma model

The torsion trace vector acts as a Proca dark photon kinetically mixed with the photon, with an effective plasma mass $m_A^2$ on the photon. Earlier runs were in the wrong physical regime (a sign error meant $\alpha_3 > 0$ was the tachyonic-instability regime, not the stable Proca regime). After correction:

| Run | $\log Z$ | $D_{\rm KL}$ | Verdict |
|---|---|---|---|
| Amplify | $-0.073 \pm 0.007$ | 0.024 nats | Null amplify |
| Suppress | $+0.66 \pm 0.05$ | 1.98 nats | Genuine suppression at decoupling corner ($m_A^2\approx 0.97, \alpha_3 \approx 0.001$) |

### References

- Barker (arXiv:2406.12826) — Barker-PGT, $\chi$ coupling
- Shapiro (arXiv:hep-th/0103093) — Shapiro derivative couplings $\zeta_{1,2,3}$
- Bahamonde et al. (2024) — propagating-torsion phenomenology programme
- An, Pospelov, Pradler (arXiv:1302.3884) — dark-photon plasma conversion

### Next steps (provisional, subject to discussion)

- **T6 parity-odd YM-PGT-CP:** tests whether parity-violation in the torsion sector opens what even-parity nonminimal couplings did not.
- **T7/T8 complete-PGT:** exhaustive quadratic enumeration. Derivations running locally.
- **Possible reorientation toward constraint-torsion + extended nonminimal couplings.** Given T4's structure vs the entire T5 family's emptiness, it may be more fruitful to broaden the *constraint-torsion + nonminimal* class — single cross-terms in the T4 mould: scalar-curvature-photon, torsion-trace-to-photon, axial-torsion-to-photon, etc. Worth your input.

---

## 3. Stability filtering in the linearised regime

To run inference on the Ricci-EM class we had to address a fundamental issue: the Padé matrix-exponential that evaluates the field equations is robust for well-conditioned systems, but the eigendecomposition-based path used previously catastrophically rejected all parameter points in the Ricci-EM prior (including physically clean ones) due to ill-conditioning. After switching to a conditioning-robust method, the approach is to compute the growth rate of the solution directly and *exclude* any parameter point where the linearised fields grow exponentially.

This choice — excluding samples with any exponential growth above a threshold — enabled all the T4 results above. But it raises a question worth discussing:

**Physical question**: In the linearised regime, an exponentially growing mode will eventually violate the linearisation assumption ($\delta g \ll g_{\rm background}$). But should all such modes be discarded? One could argue:

- Large amplification of the Gertsenshtein signal may *require* some resonant growth to accumulate — the very mechanism by which the photon channel is enhanced might look like tachyonic instability in the linearised equations.
- Some apparent instabilities may be artefacts of the linearised approximation that are cut off in the full nonlinear theory (e.g. by backreaction on the background field, or by nonlinear saturation).
- There may be a middle ground: instabilities that grow slowly enough that the linearised solution remains valid over the physical propagation length of interest, and for which the accumulated conversion is what we actually want to measure.

Currently we treat any growth rate above $\sim 0.3\,\text{s}^{-1}$ as unphysical and return $\log\mathcal{L} = -\infty$. The question is whether this threshold is too conservative, and whether some of the discarded parameter space represents genuine physics.

---

## 4. Perturbative reduction of higher-derivative theories

Some of the most physically interesting PGT Lagrangians include terms that are quadratic in the Riemann-Cartan curvature (e.g. $b_5\tilde R^2$). These generate **fourth-order** equations of motion which, naively, carry Ostrogradsky ghost modes — unphysical negative-energy degrees of freedom arising from the higher-derivative structure.

The standard resolution is *perturbative reduction*: treating $b_5$ as a small coupling, one substitutes the leading-order ($b_5{=}0$) equations of motion back into the correction terms to eliminate the higher time derivatives, leaving a second-order system.

### Methods attempted and why they failed

**(a) JLM substitution (Jaén–Llosa–Molina).**
The most direct approach: substitute the $b_5{=}0$ equations into the $b_5$-correction terms algebraically. This works when the constraint structure of the theory is unchanged at $\mathcal{O}(b_5)$. For our $b_5\tilde R^2$ PGT theory it fails because the correction *promotes constraint fields to dynamical* — in the $b_5{=}0$ theory, several torsion components are purely algebraically constrained (no time derivatives), but the $b_5$ correction adds kinetic terms for these fields. Substituting the wrong (static) equation of motion for a now-dynamical field gives incorrect correction terms and breaks the perturbative expansion.

**(b) LPS canonical analysis (Lyakhovich–Pluschchay–Sharapov) and Dirac–Bergmann.**
The principled Hamiltonian route: work out the full constraint algebra of the $b_5$ theory and classify all first- and second-class constraints. This is exact but practically infeasible: our 18-field PGT theory with a fourth-order Lagrangian generates an intractable number of constraint equations ($\mathcal{O}(10^3)$ symbolically), and the constraint structure *changes dimensionality* at $\mathcal{O}(b_5)$ — the phase space itself gains new dimensions when the promoted fields acquire dynamics.

**The key obstruction shared by both:** there is no published recipe in the literature for the *constraint-promotion* case, where a field that was non-dynamical at leading order becomes dynamical at the next order. Both approaches implicitly assume the number of dynamical degrees of freedom is fixed.

### What was implemented instead

Rather than working at the Lagrangian level, we work directly with the linearised equations of motion. The strategy:

1. **Solve the base ($b_5{=}0$) equations exactly** using the spectral solver, obtaining the leading-order field evolution $y^{(0)}(t)$.
2. **Treat the $b_5$-correction terms as a source** for the base operator. The $\mathcal{O}(b_5)$ correction satisfies $L\,y^{(1)} = S[y^{(0)}]$ where $L$ is the second-order base operator and $S$ is the correction source built from $y^{(0)}$.
3. **Solve via the Duhamel convolution integral**: the correction is $y^{(1)}(t) = \int_0^t e^{(t-\tau)A}\,S(\tau)\,y^{(0)}(\tau)\,\mathrm{d}\tau$, which reduces to a closed-form kernel evaluated from the base eigendata.

This approach is **ghost-free by construction** — the base LHS operator $L$ is always second-order, so no Ostrogradsky modes enter. It is also theory-agnostic: the same code handles any small-parameter correction without per-theory classification. Verified against exact analytic solutions (Parker–Simon FLRW to $10^{-12}$, driven-oscillator to $10^{-14}$) and against a $b_5{=}0$ reference limit.

### Open gap: validity bound for the constraint-promotion case

When $b_5$ promotes a constraint field, that field has *zero amplitude* in $y^{(0)}$ but acquires $\mathcal{O}(b_5)$ amplitude from the source. The Duhamel answer is numerically sound, but we lack a closed-form theorem bounding the validity domain of the expansion — only the heuristic $\varepsilon\cdot\omega\cdot t$ threshold from Figueras–Kovács–Yao (2025). See §Questions.

---

## 5. Practice talk

The talk covering this project has been written and is ready to schedule. Key framings:

- The Gertsenshtein effect is presented as *astrophysically useless at GR rates* before the survey question is introduced.
- The dark-photon vacuum null is framed as an *eigenmode* argument (the graviton initial condition cannot populate the mass eigenstate in which the dark photon lives) rather than a generic kinetic-mixing triviality.
- The Bayesian inference is described as mapping expected amplification across the prior, with $\log Z$ as the single summary number.

---

## Questions

### 1. Stability filtering — how conservative should we be?

We currently exclude all parameter points where linearised fields grow exponentially above a threshold. This enabled the Ricci-EM inference results. But:

- [ ] Can large amplification of the Gertsenshtein signal occur *without* some resonant growth? If the conversion mechanism is fundamentally a resonance, the very signal we are looking for might live in the excluded region.
- [ ] Are some of the apparent instabilities artefacts of the linearised approximation, stopped by nonlinear effects in the full theory?
- [ ] Is there a physically motivated threshold below which growth is acceptable — e.g. growth slow enough that amplitudes remain in the linearised regime over the propagation distance?

### 2. Perturbative-reduction validity bound for the constraint-promotion case

The Duhamel-source approach gives a numerically correct $\mathcal{O}(b_5)$ answer for $b_5\tilde R^2$ PGT, but when $b_5$ promotes constraint fields to dynamical the standard Hamiltonian methods (JLM, LPS) cannot provide the validity bound. We have only the heuristic $\varepsilon\cdot\omega\cdot t$ threshold from Figueras–Kovács–Yao (2025) and numerical cross-checks.

- [ ] Is there literature on the validity of perturbative expansions when the number of dynamical degrees of freedom changes order-by-order that you know of?
- [ ] Should we treat numerical agreement as sufficient, or is a formal bound needed before these results appear in a paper?

### 3. Campaign direction — constraint torsion vs propagating torsion

The only completed theory with posterior structure has constraint torsion + a nonminimal coupling. The entire propagating-torsion nonminimal sector (five sub-theories) is null.

- [ ] Is it worth pivoting to broaden the constraint-torsion + nonminimal class (T4-mould theories with different cross-terms) rather than continuing into propagating-torsion / parity-odd / complete-PGT?
- [ ] Or complete the parity-odd and full-quadratic enumeration first?

### 4. Initial-condition dependence in the Ricci-EM results

Initial runs give a strong suppressor (Bayes factor 0.10 for amplification). A rerun with the canonical Fourier-grid initial condition gives an apparent amplifier ($\Delta\log Z = +4.4$ nats, apparent $A_{\max} \approx 42$, but borderline perturbativity). Cross-check still in flight.

- [ ] If the verdict does flip with initial condition, what is the right choice for paper-grade reporting?

### 5. Reporting a survey of nulls

The dark-photon, Einstein-Cartan, and all five YM-PGT sub-theories are null. The nonminimal couplings (six parameters spanning the propagating-torsion sector) are all flat in posterior.

- [ ] What is the right framing for the paper — a joint Bayes factor, per-theory, or something else?

### 6. Carry-forward from 17 April

- [ ] **Scheduling the practice talk** — when and to whom? Would you like to be in the audience?
- [ ] **Will Handley** PhD application — when to email?
- [ ] **Sven Krippendorf** — outcome of any intervention?
