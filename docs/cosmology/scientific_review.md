# Scientific review before implementation (2026-09-06)

**Status:** complete. **Scope:** the eight design handoffs H1–H8 and the eleven documents
they produced, reviewed as one body of work before any implementation is delegated.
**Verdict: the programme is sound and implementation proceeds on the settled path.**

The rung order (`O0 → O1 → O2 → O4a → O3 → O4b/V`) and the milestone sequence
(`M0 → M0.5 → M1 → M2`, WS6 in parallel) are **unchanged** by this review. Two decisions
and one open question came out of it; everything else is correction and clarification.

---

## 1. Method, and what this review is not

The handoffs ran sequentially over six days, each settling part of the design. The
coherence pass (v0.52.2) reconciled them for internal consistency. This review asks the
different question: **is the physics the implementation will be built on correct, and is
it stated precisely enough to build from?**

**Physics research is trusted, not re-derived** (user decision, 2026-09-05). H3 verified
its conventions against the source papers' own benchmarks; H8 read the PSALTer sources
live; H2 worked the per-operator analysis. Re-deriving that would mostly re-verify what
prototypes already confirmed. The review's job is comprehensive orchestrator
understanding, precision of the instructions an implementer will follow, and a defensible
plan for what happens next.

Evidence base: three exhaustive document audits (physics/observables; solver + CAMB seam;
spectrum + engineering + migration), each finding quoted to `file:line`; first-hand reads
of `spectator_route.md`, `solver_design.md` §4, the dependency graph, and
`background_validity.tex`; a red-team pass over the delegation plan verified against the
repository. ~40 findings, all dispositioned in §4.

**Out of scope:** re-deriving handoff physics; the GPL/MIT release gate (#495); HPC (D4).

---

## 2. Per-rung verdicts

| Rung | Mechanism basis | Known-answer oracle | Verdict |
| --- | --- | --- | --- |
| **O0** — pass-through | No new physics | CAMB's own `C_ℓ` | **Sound.** Gate de-circularized by this review: a pass-through returns CAMB's arrays, so *machine-precision identity* is the honest test — "sub-percent" would pass by construction while hiding a unit or convention slip (`repo_reshape.md` M1a row) |
| **O1** — tabulated background | H1's audit of TorC; no TorC physics enters our package | Constant table → ΛCDM `C_ℓ` to machine precision; `TorC_rhopa.py` table → reference `H(a)`/`C_ℓ`; free `set_w_a_table` cross-check at `ϖ_r = 0.8` | **Sound.** A plumbing gate by design (R1). Risks are mechanical, and H1 §8 supplies a second independent route to the same physics at the fiducial point |
| **O2** — tensor propagation | Coupled `(h, torsion 2⁺)` block, reducing to Cembranos eq. 20 | CAMB's tensor `C_ℓ^BB` in the decoupled limit, sub-percent | **Sound, one gap.** No WS2→WS3 interface contract exists (§3.1). Sourced evolution is correctly first-class — the gate is unreachable without the anisotropic-stress source |
| **O4a** — isotropic `β` | CS/photon sector conformally invariant → `β` from a zero-mode integral; no photon solver | `class_rot` at `ᾱ = 0.1°` to `10⁻⁶`; Das et al. bound; `β = 0.277° ± 0.057°` (4.8σ) | **Sound, two conditions.** Requires `T̄ ≠ 0`, so the background-EOM residual is mandatory here and was missing (§3.2). Its cheapness rests on `n = 0` **and** `β` constant — #503 settles only the first |
| **O3** — Gertsenshtein on FRW | Eikonal amplitude system with patch averaging | Flat-space `P = sin²(κB₀D/2)` to 0.4%; Cembranos eq. 34 and `P = 8.6×10⁻¹²`; patch-averaged `𝓕` to factor ~2 | **Sound, modeling-limited.** Conventions pinned against the paper's own numbers, including two errors caught in the source. The binding uncertainty is not numerical but the assumed B-field (§5) |
| **O4b / V** | Rotation inside the LOS integral; chiral conversion and `E→V` | `class_rot`; published V-mode bounds only | **Sound, least specified.** No likelihood exists for V-modes anywhere and Cobaya's `Cl` keys cannot express them; correctly scoped as a bound |

**The architecture holds at every rung.** The spectator route is standard methodology
(five names across five subfields), its scope boundary is drawn correctly, and the one
thing the literature never does — enforce validity numerically rather than asserting it —
is the programme's own contribution.

---

## 3. What this review changed

### 3.1 One real design gap: O2 has no WS2→WS3 contract

O3's symbolic interface is specified as four numbered deliverables plus an acceptance test
(`solver_design.md` §5). **O2's is one sentence about array shapes.** There is no stated
JSON contract for the state vector `(h, h'/k, δT, δT'/k)`, no definition of the sources
`S_std`/`S_new` (which appear in the block equation and nowhere else), and no specification
of how CAMB's `π_ν` table reaches the augmented matrix — while **A8, the sourced run
against real CAMB tables, is the gate for the whole O2 stage.**

Tracked as its own issue, companion to #504. It must be settled before WS3's O2 front-end
is delegated.

### 3.2 One open physics question: does an FRW background solve the PGT field equations?

The spectator route's order-1 step — "tadpoles vanish on the background solution" — is what
licenses keeping the quadratic action whole. The order-1 coefficient **is** the background
field equation, so it vanishes only if the background solves the equations of *the theory
being expanded*. CAMB's solves **Einstein's**; ours is **PGT**. A surviving linear term is
a **source** in the perturbation EOM: the spectator modes would be driven, not free.

`background_validity.tex` proves `T̄ = 0` exact for all PGT+EM theories with non-minimal
couplings — **on flat Minkowski**, and the result it rests on finds non-trivial `T̄`
precisely for **curved** spacetime. FRW is curved. `2003.02690` gets partway: `Q = 0`
solves the pseudoscalar equation for *any* couplings, but a fully torsion-free FRW
background additionally requires `σ₃ = 0` or Einstein–Cartan.

**Handled in three parts** (user decision, 2026-09-05): declare the admissible-theory scope
as an explicit validity condition; test it per theory via the background-EOM residual
(#501), whose tolerance is *derivable* because the residual is the tadpole coefficient —
the induced source must sit far below the signal, and the expected regime is a small new
term on top of GR; and scope the FRW extension as research rather than pretending a brief
investigation settles it. **Supervisor agenda item.**

### 3.3 Two decisions

**Matrix-WKB is built alongside Magnus, not gated behind it** (user, 2026-09-05/06). The
supervisors expect WKB to be the approach that works and recommended rebuilding methods
published for analogous problems — which is what §9 does. A2b already measured the property
that matters: error *identical* at `k = 10, 100, 1000` at fixed 60 steps, where Magnus at
the same step count is useless.

**This is not a decision that Magnus loses**, and the §1 cost arithmetic is explicitly not
the evidence — it counts `expm` only, omits the assembly the same document calls dominant,
and lands at ~15 s against its own 10 s gate. Both are implemented, both are measured, and
the **bake-off decides on real numbers** with `rk-adaptive` as its mandatory measured
baseline. No candidate is discounted before a real attempt (D7). Recorded honestly: no
matrix RKWKB implementation exists anywhere, so this is a generalization of a scalar
method — the highest-risk piece of WS3, front-loaded deliberately.

**Admissible theories are scoped explicitly**, with the background-EOM residual promoted
from diagnostic to first-class per-theory gate (§3.2).

### 3.4 Supervisor intel that lowers a risk

Barker: the massless analysis was omitted from the released supplementary material **for
convenience — TorC was not interested in massless particles — not because it is hard**, and
the general algorithm is **well understood** in the literature. So "implemented numerically
nowhere" is scoped to *the released code*, and L4 becomes **find-and-implement**, not
invent: a literature search precedes implementation in the Stage-2 handoff. L4 drops in the
risk ranking accordingly.

---

## 4. Finding → disposition ledger

**A = amended at the instruction site · B = decision recorded · C = wired into a gate**

| # | Site | Finding | Cat | Resolution |
| --- | --- | --- | --- | --- |
| 1 | `solver_design.md` §5 | Acceptance test required Cembranos **eq. 24**, whose printed `/2` §4 rejected numerically (4–25× miss) — a factor 4 in `P`, inside the test meant to catch it | A | `237c33e5` — follow §4; the paper's prose/number conflict recorded |
| 2 | `solver_design.md` §4 | Two-way-coupling threshold stated at `10⁻³` and `1` in adjacent sentences | A | `237c33e5` — separated: pass / **flag** / fail |
| 3 | `solver_design.md` §4 | Seljak–Zaldarriaga cited *for* the tensor source while the cited lines call it "always negligible" | A | `237c33e5` — re-anchored to Weinberg (2004); S–Z kept as precedent for the practice |
| 4 | 7 files | `ν = 100 GHz → k ≈ 2×10²² Mpc⁻¹`, `~10²⁶` oscillations — **~10³ low** (`6.5×10²⁵`, `~10²⁹`) | A | `237c33e5` — corrected everywhere; conclusion unchanged, case strengthened |
| 5 | `primer.md` | Presented `δρ, σ` gravitational sourcing as "precisely O2" — the explicitly *unreachable* extension, in the document newcomers read first | A | `237c33e5` — O2 restated as tensor propagation; the rest marked as the road not yet taken |
| 6 | `solver_design.md` §9 | Matrix-WKB pseudocode unimplementable: `Φ` undefined where full-`Φ` is fatal; `Γ` one node vs GL4's two; GL3/GL4 mismatch | A | `237c33e5` |
| 7 | `solver_design.md` §1 | 82.6 ms splices two runs (components sum to 71) and was measured on a 30×30 block vs the design's `n≈4–10` | A | `237c33e5` — hedged; #518 to match |
| 8 | `solver_design.md` §1 | O2 budget estimate counts `expm` only, omits dominant assembly, lands ~15 s vs its own 10 s gate | A | `237c33e5` — marked a floor, not a verdict; no stepper discounted on it |
| 9 | 3 files | "No Boltzmann code has per-frequency photon propagation" — uncited, and reads as false since codes *are* `k`-resolved | A | `237c33e5` — rewritten as `k`-resolved but frequency-**integrated**; #511 evidence; falsification note; `k`-vs-`ν` added to `primer.md` |
| 10 | `observable_ladder.md` §4.3 | O3 target `P = 1.67` violates the programme's own `P ≪ 1` flag | A | `237c33e5` — marked formula-check-only, trips the flag by design |
| 11 | `solver_design.md` §4 | `modal.py:3775` "legacy precedent" cites a different augmentation | A | `237c33e5` — **removed**, not corrected; WS3 is ground-up (user) |
| 12 | `spectator_route.md` §3 | Tadpole step conditional on a check presented as a safety net; FRW case open | A+B | `bd5deccf` + new issue (§3.2) |
| 13 | `observable_ladder.md` §4.2 | O4a requires `T̄ ≠ 0` yet had no residual check; `RES` wired only into O2/O3 | A+C | `bd5deccf` — added to O4a's capability table |
| 14 | `observable_ladder.md` §4.2 | #503 gates only `n`; "`β` constant over recombination" gated by nothing | A+C | `bd5deccf` — flagged; #503 rescoped |
| 15 | `spectrum_design.md` §6.2 | "Implemented numerically nowhere" read as *invent it* | A+B | `bd5deccf` — scoped to released code; literature search now a Stage-2 deliverable |
| 16 | `solver_design.md` §9 | `5π/max\|Im λ\|` attributed to riccati; appears in neither source | A | `bd5deccf` — marked our design choice, bake-off sets the constant |
| 17 | `solver_design.md` §7 | "Memory never disqualifies — production runs land on HPC" vs **D4** | A | `bd5deccf` |
| 18 | `COSMOLOGY_PROGRAM.md` | `\|f\| ≪ 1` dimensionally ill-formed | A | `bd5deccf` — restored to `\|f\| ≪ \|F̄\|` |
| 19 | `spectrum_design.md` ×2 | `gauge_certificate` referenced as the convention to follow; `repo_reshape.md` §5.3 **drops** it | A | `bd5deccf` — renamed to the flagged-rejection convention, defined at M1a |
| 20 | `stage1_engineering_plan.md` §11 | Stale "six inline amendment notes" | A | `bd5deccf` — defers to the live ledger |
| 21 | `observable_ladder.md` §4.3 | `ν²` exclusion grounded on the paper's `n = −2` (**Faraday**) clause; our operator needs `n = +2` | A | `bd5deccf` — conclusion holds at ~5.4σ, from the other tail |
| 22 | `tests_cosmo/` | Boundary test scanned `tidalcosmo/` only; §8 requires `tests_cosmo/` too and calls it the half that matters most | A+C | `7b6f7a17` — extended, **verified by a probe that made it fail** |
| 23 | `repo_reshape.md` §8 | M0.5 fixture list asked legacy for `C_ℓ`/transfer arrays it cannot produce; omitted `validate` verdicts | A+C | `7b6f7a17` — rescoped; `C_ℓ` moves to M1a; measured scalars cut |
| 24 | `repo_reshape.md` §7 | M1's gate circular (pass-through vs sub-percent); M∥ said "exactly" where the design says sign+`10⁻⁶`, and omitted the 1 ms gate | A+C | `7b6f7a17` — M1 split into M1a/M1b; M∥ corrected |
| 25 | `repo_reshape.md` §7 | M∥ vs M0.5 "before any porting" and vs M3's ownership of `derive/` | A | `7b6f7a17` — resolved, both directions stated |
| 26 | `solver_design.md` §8/§12, #519 | WKB contingent on a bake-off | B | `237c33e5` — planned rung; bake-off decides composition |
| 27 | — | O2 has no WS2→WS3 contract | C | New issue, companion to #504 (§3.1) |
| 28 | `solver_design.md` §2.1 | CAMB `derivst` convention unverifiable in-repo (no CAMB installed; scratchpad anchor) | C | Pinned inside the A8 gate: read installed source, commit the mapping + a test. *Stale premise, 2026-09-15: CAMB 2.0.4 is installed and R-1 pinned its variable definitions from source (`conventions.md` §3); the mapping test is I-532's* |
| 29 | `solver_design.md` §9 | `[survey]`-tagged LJL bounds — the only quantitative accuracy claim for the WKB rung — unverified, no tracker | C | New tracker issue; verify on first use |
| 30 | `magnetic_field_background.md` §3 | PMF anisotropic-stress bound engaged by one paper and set aside on a disputed criticism; §1's rebuttal is an *energy-density* argument that does not answer it | C | #506 — must-adjudicate-before-O3 |
| 31 | `solver_design.md` §11 | `rk-adaptive`, the declared bake-off baseline, was never measured | C | Mandatory first row of the WS3 bake-off |
| 32 | 4 sites | Interfaces written against but never defined (§6) | C | Assigned to the handoffs that own them |

Findings 33–40 are the audits' smaller items — stale counts, `$Local*` 6-vs-8, P_max/P_final
flag semantics, `ΔN_eff` applied to possibly non-relativistic torsion — carried into the
handoff prompts that touch them rather than amended blind.

**My own false positive, recorded:** I first filed the FRW background question as a
*circularity in the argument*, then over-corrected to *settled by the thesis*. Neither was
right: the thesis settles Minkowski, FRW is open. Recorded because a mis-filed severity
wastes exactly as much of an implementer's time as a missed finding.

---

## 5. Ranked risk register

| # | Risk | Why it ranks here | Resolving test | Owner |
| --- | --- | --- | --- | --- |
| 1 | **FRW background may not solve the PGT equations** | Conditions the spectator expansion itself; open, and O4a sits outside the safe class by construction | Background-EOM residual per theory, tolerance from "induced source ≪ signal" | #501 + new scoping issue |
| 2 | **Matrix-WKB has no implementation anywhere** | Deliberately front-loaded; a generalization of a scalar method, not a port. A4's oracles (`oscode`/`riccati`) are not even installed | The bake-off, against `rk-adaptive` and A2b's k-independence | #519 |
| 3 | **PMF anisotropic-stress bound unadjudicated** | The weakest link in O3's setup: the only engaged bound is set aside on a disputed criticism, and O3's answer already spans `10⁴` on the field choice alone | Adjudicate before any published O3 bound | #506 |
| 4 | **O4a's cheapness may not survive** | Needs `n = 0` **and** constant `β`; if either fails it is O4b, and the settled rung order loses its rationale | #503 (dispersion) **and** the zero-mode EOM | #503 |
| 5 | **A8's CAMB conventions unverifiable in-repo** | The `8πG a²Σρπ` normalization sits inside O2's gate, anchored to a scratchpad copy | Read installed CAMB source; commit mapping + test | A8 / M1a |
| 6 | **`O_LL` conditioning threshold undefined** | Schur complement needs invertibility; near-degenerate cases have no stated policy | Set during Stage-2 implementation, with the accidental-degeneracy flag | #495 |
| 7 | **Spectator error never quantified** | No paper in the mixing literature bounds the test-field approximation's error; we inherit that | The four validity monitors, reported per run | WS2/WS4 |
| — | ~~L4 massless algorithm~~ | **Dropped** — Barker: standard algorithm, omitted for convenience (§3.4) | Literature search then implement | #495 |

---

## 6. Interfaces that must be defined before their consumers are built

Four cases of "module A will call this, module B will provide it" with no agreed shape —
how two sessions build incompatible halves:

1. **How a run reports that it is untrustworthy.** Five documents reference a shared
   flagged-rejection mechanism; none defines a flag's fields, severity, or how a likelihood
   consumes it — and the READMEs contradict each other on whether a flag may ever *reject*
   a sample or is purely advisory. **M1a defines it, with both severities**, because the
   spectrum module needs it as soon as it is wired in.
2. **What the solver may ask the background for** (`background/protocol.py`). **Not a gap
   to close before dispatch — it is M1a's own first deliverable**: investigating CAMB's API
   and deciding the integration is that session's work, and the protocol is its output.
3. **The Stage-1 → Stage-2 handover objects.** Field names and types unspecified; the
   unitarity-conditions status needs **four** outcomes (closed form / provably impossible /
   absent / PSALTer halted), not a boolean. Owned by the Stage-1 Python handoff.
4. **How a symbolic inequality becomes a fast per-sample test.** The massless health check
   is exported as an algebraic condition in the couplings, while the design bans symbolic
   evaluation on the per-sample path. **Genuinely unsolved** — decided in the Stage-2
   handoff rather than assumed settled.

---

## 7. Path confirmation

Unchanged: the rung order, the milestone sequence, the integration target (option iii), the
two-engine solver architecture, the two-stage spectrum architecture, and every
known-answer oracle.

Changed: WKB's sequencing (§3.3), the admissible-theory scoping condition (§3.2), and the
M1a/M1b split with de-circularized gates.

**Implementation may begin.** Wave 0 — M0 packaging (#524), M0.5 oracle freeze (#525),
PSALTer install (#526) — carries no dependency on any finding above.
