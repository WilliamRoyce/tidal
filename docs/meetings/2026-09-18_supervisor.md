# Supervisor Meeting — 18 September 2026 (DRAFT, in preparation)

**Period**: 29 August (programme pivot) to 18 September 2026.
**Status of this file**: living draft. Points are added as they arise rather than
reconstructed on the day; anything unresolved by the meeting stays here as an open item.
Moved from 11 September, when the meeting did not take place; §1.1 was added on 18 September.
The Headline and §5 still describe the period to 11 September.

---

## Headline

The cosmology programme is **design-complete, and the first implementation wave is merged**.
Eight research handoffs (H1–H8) ran between 29 August and 4 September, producing thirteen
design documents; a coherence pass reconciled them, and a scientific review before dispatch
(`docs/cosmology/scientific_review.md`) confirmed the architecture holds at every rung.

**The first wave delivered all three goals, and the third took a detour worth reporting.**
The new package is installable with a CI lane that found a real defect on its first execution;
the legacy oracle is frozen as 185 committed fixtures before any porting; PSALTer is installed,
three source questions are answered from a live install, and **its install gate now passes** —
but only after a fortnight's detour that ended in a finding about your package rather than
about ours (§3b). We also retired four legacy subcommands the new design supersedes.

**Settled since the last meeting:** the observable ladder's execution order
(`O0 → O1 → O2 → O4a → O3 → O4b/V`), the integration target (our own solver chained to
*unmodified* CAMB), the two-engine solver architecture, and the two-stage spectrum
architecture.

---

## 1. The question I most want your view on — does an FRW background solve the PGT field equations?

**Why it matters.** The whole spectator route rests on expanding the action about a
background and having the order-1 term vanish. That order-1 coefficient *is* the background
field equation, so it vanishes only if the background solves the equations of the theory
being expanded. We take CAMB's background, which solves **Einstein's** equations, while our
action is **PGT**. If a linear term survives, it is a **source** in the perturbation
equations — the spectator modes would be driven rather than freely propagating.

**What we already have.** The thesis settles the flat case:
`docs/tex/background_validity.tex` §"Background Torsion: `T̄ = 0` Is Exact" shows the Cartan
equation gives `0 = 0` for all PGT+EM theories with non-minimal couplings — **on flat
Minkowski with uniform `B₀`**. But the result it rests on (Bahamonde et al.) finds
non-trivial `T̄` precisely for **curved** spacetime, and FRW is curved.

**What we know on FRW.** From `2003.02690` (the group's own): `Q = 0` solves the
pseudoscalar torsion equation identically for *any* couplings — encouraging — but a *fully*
torsion-free FRW background additionally requires `σ₃ = 0` (k-screening) or the
Einstein–Cartan case. Outside those, theories sit in a tracking class whose effective
gravitational constant is rescaled, so CAMB's `G` is not the theory's `G`.

**Our provisional handling**, which I would like checked:

1. **Scope it.** Admissible theories are those admitting the assumed background — stated as
   an explicit validity condition of the route, with the theory judgement left to the user
   rather than silently assumed by the pipeline.
2. **Test it per theory, in-pipeline.** A background-EOM residual computed on the CAMB
   background. Because the residual *is* the tadpole coefficient, its tolerance is
   derivable rather than arbitrary: the induced source must sit far below the signal being
   computed. The expected regime is a **small new term on top of GR**, where residual and
   induced source are correspondingly small.
3. **Scope the extension as research, not resolve it now** — survey which theory classes
   provably admit a torsion-free FRW background and what settling it would require.

**Questions:** Is (1)+(2) the right posture, or does this need settling before O2? Is the
tracking class worth supporting, or is restricting to `T̄ = 0`-admitting theories the
cleaner scope? And is anyone aware of work extending the Minkowski argument to FRW?

**Where it bites first:** O4a (isotropic birefringence) *requires* a homogeneous torsion
mode `S₀(η) ≠ 0`, i.e. `T̄ ≠ 0` by construction — so it is the one rung guaranteed outside
the safe class.

### 1.1 Added 18 September — a scope decision I would like you to confirm

**The decision (provisional, 17 September).** We expand about a background with **exactly
zero torsion**. The new fields and couplings are assumed to change only the perturbations, on
the same ΛCDM background evolution — which is what "spectator" means for us. We do not need a
background torsion for the effects we are after: the operators we consider couple the
perturbations directly, and that is where parity mixing and the conversion channels come from.

**What it removes.** Any calculation with a nonzero background torsion. In particular,
isotropic birefringence from a homogeneous torsion mode (O4a above) is no longer a planned
rung — on 13 September we had already moved birefringence to parity-odd couplings among the
perturbations, which need no background of any kind. It also removes an earlier allowance for
"a small nonvanishing background torsion that stays below the tolerance".

**What it keeps.** Every perturbation-level channel, including the parity-mixing ones and the
Gertsenshtein-type conversion; and the in-pipeline background check of point 2 above, which
now asks a sharper question: does the theory allow zero torsion on the ΛCDM background? A theory
whose equations drive a background torsion is out of scope, rather than handled approximately.

**The honest tension.** This decision leans on the question at the top of this section. From
`2003.02690`, zero torsion is exact only in some classes (`σ₃ = 0`, or Einstein–Cartan); the
tracking class is excluded by construction. If the theories we care about turn out not to allow
zero background torsion on FRW, the decision has to be revisited — that is the stated condition
for reopening it.

(Separately, and not the reason for the decision: the perturbation tool we use, xPand, cannot
at present take a nonzero background for a field shaped like the torsion. We have diagnosed
why and drafted a report for its author.)

**Questions:**

1. For first results, is it right to restrict to theories that allow exactly zero background
   torsion on FRW, accepting that this excludes the tracking class and isotropic birefringence?
2. For the couplings we are most interested in — parity-odd couplings between the torsion and
   photon perturbations — do you expect a torsion-free FRW background to be allowed, or should
   we plan the extension now rather than on demand?

---

## 2. Solver direction — WKB, as you expected

Confirming that the research supports the expectation from the last meeting, and reporting
one honest caveat.

**The design does what you recommended** — find the published methods for analogous
problems and rebuild them. The WKB rung is built on Lorenz–Jahnke–Lubich adiabatic Magnus,
lifted to first-order systems with non-normal `M`, cross-checked against Ioannisian–Smirnov
closed forms, with neutrino oscillation in matter (arXiv:0803.1967) as the template and
Handley's `oscode`/`riccati` as the scalar-case prior art.

**The prototype confirms the property that matters.** On a de Sitter adiabatic band at a
fixed 60 steps, the adiabatic stepper's error is *identical* for `k = 10, 100, 1000`, while
4th-order Magnus at the same step count is useless. That k-independence is the whole point.

**Caveat: no matrix RKWKB solver exists anywhere.** `oscode` and `riccati` solve *scalar*
second-order ODEs only. So this is a **generalization of a published scalar method to
matrix systems**, not a port — higher effort and higher risk than it may sound, and
independently publishable if it works.

**Decision taken:** WKB is implemented alongside Magnus in the first solver handoff rather
than deferred behind it. Both are measured; a bake-off decides composition and handover
thresholds on real numbers, with an adaptive RK baseline as the control. No candidate is
discounted on paper estimates.

One honest caveat on evidence: the error bounds we would quote for the adiabatic method are
still `[survey]`-tagged — taken from the Lorenz–Jahnke–Lubich paper at second hand and not yet
re-derived against the primary source. They are the only quantitative accuracy claim the rung
has, so they get verified at first use rather than cited as settled (#530).

**Question:** does treating the matrix generalization as a publishable result in its own
right match how you would want it framed?

---

## 3. For Wolfgang — the massless spectrum algorithm

Thank you for the steer that the massless analysis was left out of the supplementary
material **for convenience** (TorC not being interested in massless particles) rather than
because it is hard, and that the general algorithm is well understood.

That materially changes how we plan it: our design had recorded it as "implemented
numerically nowhere", which read as *we must invent it*. It is now scoped as
**find-and-implement** — locate the published treatments, curate them, and implement the
complete algorithm — the same pattern we are using for the Schur-complement criterion.

**Question:** which references do you have in mind for the general algorithm? That would
save us a literature search and, more importantly, make sure we implement the version you
would recognize as complete.

Related, and already acted on: the released validator does **not** enforce coupling-linearity —
**measured rather than read**, since the install emits ~23 ambient messages on any theory, so
we diffed a control theory against one with a bare numeric coefficient and found the
difference empty. `NonLinearCouplings` is defined but never thrown. A fourth check,
`NonQuadraticFields`, *has* a throw site but is guarded by `ResourceFunction["PolynomialDegree"]`,
which cannot be fetched in our environment — so it is inert here too. We enforce
coupling-linearity on our side and treat the validator as load-bearing for correctness.

---

## 3b. For Wolfgang — we found it, and it is two undocumented dependencies in PSALTer

**Resolved on 11 September, by execution.** The install now passes its gate: running your
`ParticleSpectrographCTEG.m` unmodified reproduces your committed `.mx` exactly — both keys
identical, `PseudoDeterminant` back to real polynomials where it had been all zeros.

**The cause is not the engine.** We suspected a 14.2-versus-14.3 difference, because your
`.mx` header decodes to 14.2. That was wrong: we installed 14.2.1 alongside and it behaves
identically. The actual cause is that `ParticleSpectrum` depends on two Wolfram Function
Repository resources that are declared nowhere:

- `ResourceFunction["LinearlyIndependent"]` — `SymbolicNullSpace.m:31` on the master kernel,
  `IsNullVectorOfSpace.m:6` on subkernels
- `ResourceFunction["PolynomialDegree"]` — `ValidateLagrangian.m:38`, `UnresolvedPoleRow.m:11,21`

When they cannot be resolved — no network, a repository outage, or any environment where
`ResourceFunction` returns `$Failed` — **the run completes and writes its `.mx`, silently
wrong**. An unresolved call is not a Boolean, so the `If` that appends to `CommonNullVectors`
takes neither branch; every `SymbolicNullSpace` returns `{}`; no gauge symmetry is identified;
the compensator is zero; `Det[…]` is identically zero for gauge-singular blocks; and
`UnmakeSymbolic.m:75` divides by it. That is the first `Power::infy`, and every pseudo-
determinant downstream is zero. `NonQuadraticFields` never fires either — a cubic Lagrangian
is accepted.

**We think this is a deterministic trigger for your "Known bugs" item 1** — *"a sporadic error
where some of the gauge symmetries are not identified"*. An environment that cannot reach the
repository loses gauge identification every time rather than occasionally, through exactly the
code path that item describes.

**Why it is invisible.** `PSALTer.m:18-20` redefines `Message` to `Null` on subkernels, so the
`ResourceObject::notfname` from the subkernel call site is never seen. On the master only three
appear, and the run carries on.

**What we did about it, and what we did not.** We registered the genuine definitions locally —
`LinearlyIndependent` from the `ResourceFunctionHelpers` paclet Wolfram itself ships, and
`PolynomialDegree` from its repository notebook, hash-pinned. **No edit to PSALTer**, no
workaround inside your code, and the pinned revision is untouched. With them present the same
inputs give Booleans at both call sites, CTEG's source constraints come back carrying the
formulation's 21 generators, and the spectrograph renders cleanly.

**We have a draft issue written for your repository and have not filed it** — it is yours to
look at first: `docs/cosmology/psalter_543_upstream_issue.md`, with the reproduction, the
standalone shape of the failure, and a suggested fix (declare the two resources, or fail loudly
at load rather than silently at runtime).

**Questions:** are those two resources meant to be hard dependencies? And would you like the
issue filed as written, adjusted, or not at all — it is your package, and we would rather you
saw it before anyone else did.

**One honest note on how we got here.** We first reported the resources as "ruled out by test",
because an earlier session substituted both and got an identical verdict. Re-reading that
substitute showed it reproduced the disabled behavior rather than restoring it — our stand-in
returned a list where the real function returns a scalar, so the comparison it fed never
evaluated and the guard stayed off either way. The control had not discriminated. That is worth
saying because it is the kind of mistake that closes a question prematurely, and it cost us two
days of chasing the engine.

**Also worth knowing:** `Method` is inert on v2.0.2 — zero `OptionValue@Method` sites against
five for `MaxLaurentDepth`, and `"Easy"`, `"Hard"` and a deliberately invalid value all return
byte-identical results with identical timings. We report `ParticleSpectrum` wall time without a
Method qualifier as a result.

## 4. The weakest link in the Gertsenshtein rung — the primordial magnetic field

O3 cannot be posed without an assumed background B-field, since the mixing is *linear* in
it. Two things concern me:

- **The assumption is worth ~10⁴ in the answer.** Published choices run from 47 pG to 5 nG,
  and `P ∝ B₀²`. Any bound we quote must carry its assumed field.
- **The backreaction question is unadjudicated.** Surveying the conversion literature, *not
  one* paper justifies neglecting the field's effect on the expansion history — each takes
  a fixed classical background and imports an observational upper bound. The one paper that
  engages the relevant bound (Caprini–Durrer anisotropic-stress limits) sets it aside on
  the strength of a published criticism. Our own justification is an **energy-density**
  argument (`r_B ≈ 10⁻⁷ B₋₉²`, hence `ΔN_eff ≲ 10⁻⁵`), which does not answer an
  *anisotropic-stress* bound.

**Question:** is the energy-density argument sufficient for our purposes, or should we
adjudicate the anisotropic-stress bound before quoting any O3 result? Currently flagged as
must-resolve-before-publication.

---

## 5. Status, briefly

- **Design documents:** thirteen, under `docs/cosmology/`, with `docs/COSMOLOGY_PROGRAM.md`
  as the operational record (decisions register, ladder, workstreams, wave board).
- **Package:** `tidalcosmo/` is now a real installable package beside legacy `tidal/` — two
  console scripts, `camb`/`cobaya` extras, a CI lane on the integration branch. New code
  never imports legacy, test-enforced. (That guard inverts into a blocker at the final
  rename, which we found by asking what future change makes each guard wrong; it is
  scheduled for deletion rather than adaptation.)
- **First implementation wave — merged, one goal unmet:** packaging ✅; the legacy oracle
  frozen as 185 fixtures before any porting ✅; PSALTer installed and its three live-source
  questions answered ✅; **its Tier-1 install gate passes as of 11 September** (§3b) — the
  install is certified on Wolfram 14.3.0 × PSALTer `bb45adb0` × local registration of the two
  Function Repository resources. The second wave is planned from here.
- **Approach to delegation:** self-contained handoff prompts to separate sessions, each with
  quantitative success criteria stated before code, merged centrally against a checklist.
  Working well; the two things that bit us were a gate nobody could run and a verification
  step promised in a prompt but never copied onto the checklist that gets executed.

---

## Open items carried in

- The `ν⁰` vs `ν²` frequency scaling is **per-operator**, not a single number. Only `n = 0`
  operators can *explain* the 4.8σ birefringence signal; `ν²` operators can only be
  bounded. Derivation in progress.
- The Chern–Simons couplings need **bare `A_μ`** handling, which the pipeline does not yet
  support — an unsolved problem rather than a configuration step.
- Licensing for code derived from the PSALTer/supplementary sources: a release gate, not an
  implementation one. Attribution to be settled at publication. Concretely, we have now
  committed two small `.wxf` outputs from your supplementary materials (692 B and 2,874 B,
  with provenance and hashes) so our reader can be tested without Wolfram installed. Our
  position is that these are *data outputs of a computation* rather than "the Program" — but
  it is a position, not a settled conclusion, and we would rather hear your view now than at
  publication.
- **A positive result on cost, since it answers our own go/no-go:** CTEG — your 21-generator
  PGT — completes in **407 s** on our install, with ~72 % of that in field declaration and
  decomposition. So adding couplings to an existing field content is far cheaper than adding
  fields, and the "minutes, not hours" regime the design assumed holds.
