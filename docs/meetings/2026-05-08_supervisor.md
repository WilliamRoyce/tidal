# Supervisor Meeting — 8 May 2026

**Period**: 17 April (last meeting) to 8 May 2026

---

## Summary

Four workstreams since the last meeting:

1. **T4 Ricci-EM is the live wire.** First completed theory with posterior structure. Under the v5 IC it is a strong suppressor (Bayes factor 0.10 against amplification); under the canonical IC the rerun appears to amplify ($\Delta \log Z = +4.4$ nats). Verdict held until the suppress cross-check lands.
2. **Survey progress and roadmap.** Stage A v5 (T1 dark-photon-plasma), Stage B (T2 Einstein-Cartan), Stage D1 (T4 Ricci-EM), Stage D2.0–D2.3 (the entire propagating-torsion nonminimal sector: Bahamonde / Barker / Shapiro / full T5 9-D) are complete. T6 (parity-odd), T7/T8 (complete-PGT), and EH (higher-curvature EM) remain.
3. **Methodology updates that changed what's accessible.** Stability-guard refactor (T4-class theories now usable); perturbative-reduction v6 (constraint-promotion theories like $b_5\tilde R^2$ PGT now usable, ghost-free by construction); CDT sign-convention flip (Stage A now in the physical Proca regime); unified $\log_{10} A$ Bayes-factor framework.
4. **Practice talk drafted.** Script finalised; framing tightened by audience-test of the narrative.

---

## 1. T4 Ricci-EM (the live wire)

**Theory**: `examples/torsion_gertsenshtein_nonminimal/theory.toml`

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}R + \alpha_i\,I_i + \delta_1\,\tilde R_{\mu\nu}\,F^{\mu\nu} - \frac{1}{4}F_{\mu\nu}F^{\mu\nu}$$

The three torsion-mass invariants $\alpha_{1,2,3}$ plus a single nonminimal coupling $\delta_1$ between the Ricci-Cartan tensor and the photon field strength.

**Crucially, T4 has *constraint* (non-propagating) torsion** — there is no kinetic term for torsion, so it carries no dynamical modes and acts purely as an auxiliary field. $\delta_1$ is the only term that *connects* torsion to the Gertsenshtein channel; without it, torsion is structurally decoupled (the takeaway from the 17 April meeting). This makes T4 the only completed theory in the *constraint torsion + nonminimal coupling* quadrant of the survey — see §2.

### Stage D1 v5 (paired runs, $k_{\rm IC}=2.0$)

| Run | Job | $\log Z$ | joint $D_{\rm KL}$ | Posterior signal |
|---|---|---|---|---|
| Amplify | 28520217 | $-2.26 \pm 0.07$ | 1.79 nats | $A_{\max}=1.26$; Bayes factor 0.10 vs null (model **disfavoured 10:1** for amplification) |
| Suppress | 28519675 | $+15.92 \pm 0.13$ | 8.91 nats | $A_{\min} \approx 4\times 10^{-12}$ at MAP; valley reaches $\sim 5\times 10^{-9}$ across $(\alpha_i, \delta_1)$ |

- $\delta_1$'s marginal $D_{\rm KL}$ is **96% of the joint** in amplify — a single coupling carries the entire structure. In suppress, $\delta_1$ leads (0.87 nats) but $\alpha_{1,2,3}$ each contribute $\sim 0.25$–$0.30$ nats: the deepest suppression valley needs all four parameters to coordinate.
- Suppression depth is far below the analytic estimate ($\sim 10^{-3}$). The valley signature is consistent with destructive interference, but the mechanism is not yet pinned down.
- **This is the first non-trivial result of the survey** — every other completed theory has been null.

### Phase 6.C canonical-IC rerun (in flight, deferred)

Same Lagrangian, but with $k_{\rm IC} = 2\pi/L$ (the canonical Fourier-grid choice, $L=100$) and the Hwang–Noh stability gate.

- **Amplify rerun (28789579, INTR, completed 2026-05-03)**: $\log Z = +2.135 \pm 0.059$, a $+4.4$ nat shift from v5. New MAP at $\delta_1 \approx +1.94$, apparent $A_{\max} \approx 42$.
- **Caveat**: $\sim$26% of the posterior has $P_{\max} > 0.3$, so perturbative validity is borderline.
- **Suppress cross-check (28799598, INTR)**: still running. **Verdict held until it lands.**

If the conclusion really does flip with IC choice, that itself is a methodological finding — and a question for the meeting (see §Questions).

### Literature

| Claim | Citation | Reference |
|---|---|---|
| Ricci $\times F$ as the sole single-coupling channel | This work — informed by 17 April meeting | T4 derivation |
| Constraint vs propagating torsion in PGT | Hehl et al. 1995 | _Phys. Rep._ 258, 1 |
| Padé matrix-exponential for ill-conditioned linear systems | Higham 2008 | _Functions of Matrices_ |

---

## 2. Survey progress and roadmap

### Theory partition

The linearised landscape splits naturally on two axes — *whether torsion propagates* and *whether nonminimal couplings are present*:

| | constraint torsion | propagating torsion |
|---|---|---|
| **minimal (no nonminimal coupling)** | T2 Einstein-Cartan — null | T3 $R^2$-PGT — $b_5$ decouples structurally |
| **nonminimal coupling present** | **T4 Ricci-EM — structure** | T5 YM-PGT family (Bahamonde, Barker, Shapiro, full 9-D) — all null |

Plus the phenomenological/effective class (T1 dark-photon-plasma, plasma Gertsenshtein) and the not-yet-completed parity-odd / complete-PGT / higher-curvature-EM theories (T6, T7, T8, EH).

Within T5 we deliberately tested **nested sub-theories** prominent in the literature before the full 9-D run: Bahamonde $\subset$ Barker $\subset$ Shapiro $\subset$ full T5. Each layer adds one or more nonminimal couplings; each layer was null; the full 9-D confirms there is nothing hiding in the higher-dimensional joint.

**Emerging pattern (worth raising explicitly):** the only quadrant with structure is *constraint torsion + nonminimal coupling*. Propagating torsion is inert across the entire nonminimal subspace tested (couplings $\delta_1, \chi, \zeta_{1,2,3}$, plus the kinetic $\xi$). The kinetic structure that gives torsion its own modes appears to dilute or close the cross-channel.

### Roadmap

| Stage | Theory class | Free dim | Status | Verdict |
|---|---|---|---|---|
| A v5 | T1 dark-photon-plasma (effective) | 4 | Done | Null amplify; real suppression at decoupling corner |
| B | T2 Einstein-Cartan (minimal) | 3 | Done | Null — torsion structurally decouples |
| (deferred) | T3 $R^2$-PGT ($b_5$ minimal-quadratic) | 4 | Deferred | $b_5$ decouples from TT channel structurally |
| D1 v5 | T4 Ricci-EM (nonminimal, constraint torsion) | 4 | Done | **Strong suppressor; structure on $\delta_1$** |
| D1 6.C | T4 canonical-IC rerun | 4 | **In flight** | Held — apparent amplifier; awaiting suppress cross-check |
| D2.0 | T5-Bahamonde ($\beta, \xi, \delta_1$) | 5 | Done | Null |
| D2.1 | T5-Barker ($+\chi$) | 6 | Done | Null — $\chi$ inert |
| D2.2 | T5-Shapiro ($+\zeta_{1\text{–}3}$) | 8 | Done | Null — $\zeta_i$ inert |
| D2.3 | T5 full (9-D) | 9 | Done | Null — all 6 nonminimal couplings inert |
| **D3** | T6 YM-PGT-CP (parity-odd) | $\sim$22 | **Next** | Pending derivation; HPC submission planned |
| **E1** | T7 Complete-Even-PGT | $\sim$20 | Pending | Wolfram derivation in progress |
| **E2** | T8 Complete-Odd-PGT | $\sim$30+ | Pending | Wolfram derivation in progress |
| EH | Einstein–Maxwell + higher-curvature EM | TBD | Pending | Blocked on xAct $(F\!\cdot\!F)^2$ issue |

### D2 sub-stage results (compact)

| Stage | Theory | Amp $\log Z$ | Sup $\log Z$ | Bayes factor | Inert couplings |
|---|---|---|---|---|---|
| D2.0 | T5-Bahamonde | $+0.616$ | $-0.449$ | 2.90 | $\delta_1$ |
| D2.1 | T5-Barker | $+0.618$ | $-0.447$ | 2.90 | $\delta_1, \chi$ |
| D2.2 | T5-Shapiro | $+0.612$ | $-0.615$ | 3.41 | $\delta_1, \chi, \zeta_{1\text{–}3}$ |
| D2.3 | T5 full | $+0.6150$ | $-0.6146$ | 3.42 | $\delta_1, \chi, \zeta_{1\text{–}3}$ |

**Key claim**: the entire propagating-torsion nonminimal sector (Barker $\chi$, Shapiro $\zeta_i$, Bahamonde $\delta_1$ in T5 form) is **Gertsenshtein-neutral** — these couplings do not open conversion channels. Posterior shape is dominated entirely by the $\beta_{1\text{–}3}, \xi$ stability boundary; the *interesting* couplings are flat (max marginal $D_{\rm KL} < 0.06$ nats in D2.3 across all six nonminimal parameters).

### Stage A v5 (T1 dark-photon-plasma, post sign-convention flip)

| Run | Job | $\log Z$ | $D_{\rm KL}$ | Verdict |
|---|---|---|---|---|
| Amplify | 28474676 | $-0.073 \pm 0.007$ | 0.024 nats | Null amplify; previous "informative" v4 amplify was ghost contamination |
| Suppress | 28477675 | $+0.66 \pm 0.05$ | 1.98 nats | Real suppression at decoupling corner ($m_A^2\approx 0.97, \alpha_3 \approx 0.001$) |

The v4–v5 consistency of the suppress $D_{\rm KL}$ confirms the decoupling corner is a genuine physical mechanism, not a ghost.

### References

- Barker (arXiv:2406.12826) — Barker-PGT, $\chi$ scalar
- Shapiro (arXiv:hep-th/0103093) — Shapiro derivative couplings $\zeta_{1,2,3}$
- Bahamonde et al. (2024) — propagating-torsion phenomenology programme
- Holdom 1986, _Phys. Lett. B_ 166, 196 — kinetic mixing benchmark
- An, Pospelov, Pradler 2013 (arXiv:1302.3884) — dark-photon plasma conversion

### Next steps (provisional, subject to discussion)

- **D3** (parity-odd YM-PGT-CP): direct extension of the T5 null. Tests whether parity-violation in the torsion sector opens what nonminimal even-parity couplings did not.
- **E1/E2** (complete-even / complete-odd PGT): exhaustive quadratic enumeration. After D3 either gives null or signal, E1/E2 settle the quadratic question.
- **EH**: higher-curvature EM (Euler–Heisenberg style $F^4$). Blocked on a Wolfram-side xAct issue with $(F\!\cdot\!F)^2$ distribution; separate fix.
- **Possible reorientation toward more constraint-torsion theories.** Given T4's structure vs T5's emptiness, it may be more fruitful to broaden the *constraint-torsion + nonminimal-coupling* class — single cross-terms in the T4 mould, e.g. scalar-curvature-photon, torsion-trace-to-photon, axial-torsion-to-photon — than to keep pushing depth into propagating-torsion. Worth your input.

---

## 3. Methodology updates that changed what's accessible

Three changes mattered for results — what physical class each unlocked:

**(a) Stability-guard refactor — T4-class theories now accessible.**
The old guard rejected $\sim$100% of T4 prior samples (including physically clean points) because the eigendecomposition path failed at high condition number. Replaced with a direct matrix-exponential growth probe (Padé scaling-and-squaring), robust irrespective of conditioning. T4 inference was *impossible* before; D1 v5 is the first run that produced meaningful posterior structure on this class. Growth rates above threshold are tagged tachyonic and excluded from the likelihood — physically justified because the perturbative linearisation breaks once amplitudes reach the background scale. Per-sample rejection metadata flows through to the corner plots so the unstable region is visible, not invisible.

**(b) Perturbative-reduction v6 — constraint-promotion theories accessible.**
Theories like $b_5\tilde R^2$ PGT have constraint fields ($h_4, h_7, h_9$) that get *promoted to dynamical* by the small-parameter correction. The earlier Ostrogradsky-style reduction created singular mass matrices and ghost artefacts. The v6 idea: solve the unmodified base equations, then treat the small-parameter correction as a *source* in a Duhamel integral against the closed-form kernel of the base operator. Ghost-free by construction; theory-agnostic; no symbolic post-processing. Validity is monitored by $\varepsilon \cdot \omega \cdot t$ following Figueras–Kovács–Yao 2025. Verified against Parker–Simon FLRW to $10^{-12}$ and against driven-oscillator analytics to $10^{-14}$.

There is, however, a research gap on the validity-bound for the constraint-promotion case specifically — see §Questions.

**(c) CDT sign-convention flip — Stage A in the physical regime.**
The torsion-trace mass invariant was previously entered with the wrong overall sign: $\alpha_3 > 0$ corresponded to the tachyonic-spatial-trace regime, not the stable Proca dark-photon regime. Fixed; Stage A v5 reruns (28474676, 28477675) are in the physically intended regime. Older Stage A runs are archived as valid simulations of the wrong physics regime.

**(d) Unified $\log_{10}A$ Bayes-factor framework.**
Amplify and suppress runs now share a single derived metric $A = P_{\max}/P_{\rm GR}$, reported in $\log_{10}$. This is what enables paired-Bayes-factor comparison across the survey table. A previously-silent sign bug in the minimize likelihood briefly inverted the metric, making real suppressions look like nulls; now fixed, with the Stage A v5 suppression as the first correctly-interpreted non-trivial result.

---

## 4. Practice talk

Script finalised; the audience-test of the narrative tightened a few framings worth flagging:

- Lead with *why* the Gertsenshtein effect is astrophysically useless before the survey question, rather than after.
- The dark-photon vacuum null lands better as an *eigenmode* statement (the IC's mass eigenstate is orthogonal to the dark-photon direction) than as a generic "Holdom triviality" statement.
- The inference layer story works when $\log Z$ is described as "expected amplification under the prior", not as a generic Bayesian quantity.

---

## Questions

### 1. Validity bound for perturbative-reduction v6 in the constraint-promotion case (anchor)

When the small-parameter correction promotes a base-theory *constraint field* to dynamical (e.g. $b_5\tilde R^2$ PGT, where $h_4, h_7, h_9$ become dynamical at $\mathcal O(\varepsilon)$), neither LPS (Lyakhovich–Pluschchay–Sharapov) nor JLM (Jaén–Llosa–Molina) gives a clean reduction — both throw on our worked example. The v6 Duhamel-source approach gives a numerically correct $\mathcal O(\varepsilon^1)$ answer, but we don't have a *theorem* bounding the validity domain — only the heuristic $\varepsilon\cdot\omega\cdot t$ threshold from Figueras–Kovács–Yao 2025.

- [ ] Is it worth a focused literature search / a fresh derivation to nail this down before claiming v6 covers the class?
- [ ] Or do we treat the numerical agreement (Parker–Simon to $10^{-12}$, driven-oscillator analytics to $10^{-14}$, $b_5{=}0$ baseline regression) as the validity test in itself?

### 2. Constraint torsion vs propagating torsion — should the campaign reorient?

The only theory with structure has *constraint* torsion + a single nonminimal coupling (T4). The entire propagating-torsion nonminimal sector (T5 Bahamonde / Barker / Shapiro / full 9-D) is null.

- [ ] Is it worth pivoting HPC budget toward a class of *constraint-torsion + extended nonminimal couplings* — multiple T4-mould theories with different cross-terms — rather than continuing depth into propagating-torsion (T6, T7, T8)?
- [ ] Or run T6 / E1 / E2 first to be sure parity-odd / complete-PGT doesn't change the verdict?

### 3. Phase 6.C IC-dependence (preview, not headline)

v5 IC gives T4 as a strong *suppressor* (Bayes factor 0.10 for amplification); canonical IC gives an apparent *amplifier* ($\Delta\log Z = +4.4$ nats). Verdict held until the suppress cross-check completes.

- [ ] If the conclusion really does flip with IC choice, what's the right basis for paper-grade Bayes-factor reporting?

### 4. Reporting a survey of nulls

Stages B and D2.0–D2.3 are all null. The stability-boundary structure dominates everything visible; nonminimal couplings ($\chi, \zeta_{1\text{–}3}, \delta_1$ in the T5 form) are flat in posterior.

- [ ] For the paper, do we report a joint Bayes factor across the nested survey, or per-theory?
- [ ] Is "$B \approx 1$ across 5 nested sub-theories of propagating-torsion YM-PGT" the right framing of *strong evidence of no Gertsenshtein channel from the propagating-torsion nonminimal sector*?

### 5. Carry-forward from 17 April

- [ ] **Will Handley** PhD application — when to email?
- [ ] **Sven Krippendorf** — outcome of any intervention?
- [ ] **Practice talk** — anything outstanding on emphasis now that the script is finalised?
