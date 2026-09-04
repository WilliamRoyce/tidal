# Stage-1 engineering plan — PSALTer install, spectrum branch, exporter (H8 / WS6)

**Status:** study complete 2026-09-04 (H8). No software installed, no Wolfram run, no
pipeline code written — this document is the recommended route for the session that does.
**Design authority:** `docs/cosmology/spectrum_design.md` (H6). This document does not
re-open its decisions; it engineers them, and flags the six places where a live reading of
the PSALTer sources contradicts an assumption H6 made from the papers alone (§0).
**Tracking:** #495 (WS6), #360 (superseded PSALTer tracker), umbrella #488, decision D6.
**Findings filed:** #521 (`Method` inert), #522 (validator coverage), #523 (private-symbol
harvest surface) — each records the action it implies for the implementing session.
**Design amended 2026-09-04:** `spectrum_design.md` carries inline amendment notes at §4.1,
§4.5, §6.1 and §14, pinned to PSALTer v2.0.2 @ `bb45adb0` (§11 decision 1).
**Evidence:** PSALTer v2.0.2 at commit `bb45adb0`, and the four
`wevbarker/SupplementalMaterials-*` companions, read on 2026-09-03/04. Re-fetch with
`scripts/research/psalter_stage1/fetch_reference_sources.sh` (pins the same revisions).

## TL;DR

**Problem.** H6 settled *what* Stage 1 must hand to Stage 2 (the §6.1 contract). Nobody had
worked out *how*: PSALTer was not installed, the install had never been verified against a
published answer, no code turns a TIDAL theory into PSALTer input, no exporter emits the
contract, and the one number that decides whether the two-stage architecture survives — the
symbolic cost on a PGT+EM-sized theory — had never been measured.

**Design.** Four pieces, in dependency order: (1) a committed install script plus a
three-tier verification protocol whose cheapest tier replays Barker's own published input
and diffs against his own committed result; (2) a spectrum derivation branch that
linearizes **inside** the PSALTer session using substitution rules — the author's own
documented pattern, proven on PGT by his files — rather than the two-session route H6
assumed; (3) our own exporter, which is genuinely new code because PSALTer writes no WXF
and populates only two association keys; (4) a cost measurement whose framing must be
corrected before it is run, because `Method` is a **dead option** on the pinned release.

**Next step.** Orchestrator schedules the implementing session. It is a single sequential
Wolfram workstream (one license) with substantial Python work available in parallel; §10
gives the step order and commit boundaries.

**Done when** the §7 gates pass: a published spectrograph reproduced including its massless
condition, a Stage-1 export emitted for a campaign theory with explicit `J^P` labels, the
symbolic cost measured and reported with a go/no-go verdict, the install reproducible from
a committed script, and the two small `.wxf` files wired as unit-test fixtures.

---

## 0. Six findings that change H6's instructions

These come from reading the pinned PSALTer sources and the supplemental repositories, not
from the papers. Each one alters what the implementing session should do.

### 0.1 `Method→"Hard"` is a dead option on the pinned release (#521)

`Options@ParticleSpectrum` (`Sources/ParticleSpectrum.m:116`) declares
`Method->"Easy"`, and **`OptionValue@Method` is never read anywhere in the package** — the
only other occurrences of the symbol are the identical `Options@ParticleSpectrumActual` at
line 22 and an error message that is never thrown (§0.2). H6 §4.5 instructs Stage 1 to run
`ParticleSpectrum[…, Method→"Hard", MaxLaurentDepth→1]`, and H8's Task 4 is phrased as
measuring `Method→"Hard"` cost.

**Consequence.** On PSALTer v2.0.2 the option is inert: passing it changes nothing, and the
measurement is simply "the cost of `ParticleSpectrum` on a PGT+EM theory". Three actions
for the implementing session, in this order:

1. Pass `Method→"Hard"` anyway (harmless, and correct if a later release wires it up) but
   **do not report the number as an Easy-vs-Hard comparison**.
2. Confirm the inertness on the live install — run one small theory both ways and check the
   results and timings are identical. A static read can miss a `Method` consumed through
   `xAct`xPlain`` or an option-inheritance path.
3. Report the measurement as "`ParticleSpectrum` wall time, PSALTer v2.0.2, `Method` inert",
   never as an Easy-vs-Hard figure. The discrepancy against the design's phrasing is filed
   (#521) and **the design is amended at the point of instruction** — `spectrum_design.md`
   §4.5 and §14.1, pinned to `bb45adb0` (§11 decision 1). If the WIP `psalter.tar.gz`
   snapshot or a newer upstream commit implements `Method`, that is a finding worth carrying
   back to the author (D6 relationship), and the pinning is what keeps the amendment honest
   when it happens.

`MaxLaurentDepth` by contrast **is** consumed (validated at `ParticleSpectrum.m:29`, passed
into `ConstructMassiveAnalysis`/`ConstructMasslessAnalysis`), as is `MasslessSpectrum`,
which gates the entire massless analysis (`ConstructMasslessAnalysis.m:16`) — leave it at
its `True` default, since H6 §6.2 makes the massless sector load-bearing for us.

### 0.2 PSALTer does not enforce the coupling-per-operator contract — we do (#522)

`ValidateLagrangian.m` defines six messages but throws only four:

| message | thrown? | fires on |
| --- | --- | --- |
| `Zero` | yes | `PossibleZeroQ` on the Lagrangian |
| `UnknownCoupling` | yes | a symbol that is not a `ConstantSymbolQ` and not an `xTensorQ` |
| `UnknownField` | yes | a tensor with no `xAct`PSALTer`<name>`` context — i.e. not `DefField`ed |
| `NonQuadraticFields` | yes | `PolynomialDegree > 2` in the fields |
| `NonLinearCouplings` | **never** | — (message defined at line 6, no throw site in the package) |
| `ParityOdd` | **never** | — |

So a bare numeric coefficient (`−¼F²` with no declared coupling) **passes silently**:
numbers are not `Variables`, so nothing rejects them. H6 §4.1 says PSALTer "validates —
`ValidateLagrangian.m`, `EnsureLinearInCouplings.m`"; the second file exists but is an
internal null-vector helper under
`ConstructSourceConstraints/ConjectureNullSpace/` that clears denominators
(`Together`/`Denominator`/`FullSimplify`), not a Lagrangian validator.

**Consequence.** H6's "reject, never auto-assign" rule is entirely **our** validator's job,
and the design's expectation that a live install would demonstrate the failure is answered:
it will not. The live probe still has value — run it to confirm silence, and to capture the
exact `UnknownCoupling` text (which *is* thrown, and is the closest analogue) for our own
`error_with_hint` wording. This makes the reject rule *more* important, not less: without
it, a user's un-parameterized term is silently absorbed and every downstream verdict is
about a theory the user never declared.

### 0.3 The TorC input is published verbatim — and so is a result oracle

`wevbarker/SupplementalMaterials-2506b` is the TorC companion repository.
`ParticleSpectroscopy/system-tests-Qtorsion/ParticleSpectrographCTEG.m` (2.7 KB) is the
literal `ParticleSpectrum[…]` call behind the paper, in the **tetrad (PGT) formulation** —
`DefField[SpinConnection[-a,-b,-c], Antisymmetric[{-a,-b}], …]` plus
`DefField[TetradPerturbation[-a,-b], …]`, `TheoryName→"CTEG"`, `MaxLaurentDepth→1`,
`AspectRatio→Portrait`, `ShowPropagator→False` — i.e. it reproduces
`ParticleSpectrographTorCPGT`, the 21-generator formulation.

Two details worth their weight:

- It declares **`MuLambda` as a single constant symbol** (alongside `Mu`, `Nu`, `Lambda`,
  `MPlanck2`). That is Barker's own handling of the `μλ` product in `eq:TorCLagrangian` —
  independent confirmation that the per-operator coupling dictionary H6 §4.1 requires is
  the normal way to write these theories, not a workaround we invented.
- The repository also commits `ParticleSpectrographCTEG.mx` — a `DumpSave` of the result
  association. **This is a direct oracle**: load it and diff key-by-key against our own
  run's association, which is a far stronger check than reading a PDF.

No ECT-formulation input is published. Since H6 §4.4 makes ECT (post-Riemannian, fields
`h/T/a`, 15 generators) our primary oracle, that script is ours to author — modeled on
2506b's `PoincareGaugeTheory.m` + `Linearise.m`.

### 0.4 PSALTer writes no WXF, and populates only two association keys (#523)

`ParticleSpectrum.m:87-101` writes exactly `WaveOperator` and `PseudoDeterminant` into the
association named by `TheoryName`, via
`UpdateTheoryAssociation[…]`, which also `DumpSave`s the whole association to
`ParticleSpectrograph<Name>.mx` in `$WorkingDirectory` on every update
(`UpdateTheoryAssociation.m:14`). The spectrograph PDF is written separately through
`UsingFrontEnd@Export[…]`. **There is no WXF writer in the package.**

The richer keys H6 decoded from the committed `.wxf` files —
`UnitarityConditions`, `ComponentSourceConstraints`, `ComponentSaturatedPropagator` — must
therefore come from an **uncommitted curation step** in the polology work, which explains
H6's observation that key names vary between exports (they were assembled by hand, not by a
single writer). The genuine harvest surface inside the kernel is the pair of association
keys plus the private globals that `ParticleSpectrum` fills en route:
`$LocalSummaryOfTheory`, `$LocalSourceConstraints`, `$LocalWaveOperator`,
`$LocalPropagator`, `$LocalMasslessSpectrum`, `$LocalSpectrum`, `$LocalUnresolvedPoles`,
`$LocalOverallUnitarity` (all in `xAct`PSALTer`Private``; see `ParticleSpectrum.m:74-81`).

**Consequence.** H6 §13's "our own Wolfram-side exporter … exists nowhere and must be
written" is confirmed, and its input is now precisely identified. Two engineering
consequences: the exporter reads private globals, so it must be pinned to a PSALTer commit
and guarded by a fixture test that fails loudly if a name moves; and the `.mx` dump gives
us a zero-cost second serialization of every run for free.

### 0.5 Both linearization architectures are author-documented; the single-session one is proven

The plan going into this study assumed a two-session route (xAct session does the
post-Riemannian expansion, PSALTer session does the spectroscopy) because PSALTer owns its
own flat `M4`/`G`/`CD` with reserved indices `a`–`z` and forbids tampering. The README
documents **both** routes explicitly, and Barker's own PGT work uses the first:

- **Single session, substitution rules** (`SupplementalMaterials-2506b`
  `PoincareGaugeTheory.m` + `Linearise.m`): define the nonlinear objects — tetrad `H`,
  inverse tetrad `BField`, Riemann–Cartan `R`, torsion `T` — as ordinary `DefTensor`s on
  PSALTer's own `M4`, give them `MakeRule` definitions in terms of the perturbation fields,
  then linearize by multiplying the fields by a bookkeeping constant
  (`PerturbativeParameter`), `Series`-expanding to second order, and canonicalizing with
  `xAct`PSALTer`Private`ToNewCanonical`. The README's own words: "one can always obtain the
  correct linearisation by reinterpreting the model as a field theory on flat spacetime:
  this approach is guaranteed to work because particle physics is not actually sensitive to
  geometry per se."
- **Two sessions, copy the expression across**: also documented — "prepare the linearised
  expression in a separate *xAct* session, taking advantage of the full geometric
  interpretation, and then copy the resulting *Wolfram Language* expression directly into
  the PSALTer session. Care must be taken in this case to ensure that the indices are
  correctly matched."

**Recommendation: single-session as primary**, two-session as the documented fallback. It
removes an interchange format, an index-matching hazard, and a whole kernel invocation; it
is the route with a working PGT precedent we can read; and `Series` in a bookkeeping
parameter is exactly the ε-promotion idiom this project already uses elsewhere.

### 0.6 Inkscape is effectively unused

`$InkscapePath="inkscape"` survives in `PSALTer.m:128`, but every `Vectorize` call site is
commented out (`ConstructSpectrograph.m:49-51`, `FieldMosaic.m:10`), and `Vectorize.m`'s
own Inkscape invocation is commented out too. The spectrograph PDF goes out through
`UsingFrontEnd@Export[…]`.

**Consequence.** Install Inkscape opportunistically (one apt line, warn-and-continue on
failure) but do not treat it as a prerequisite. The real headless question is whether
`UsingFrontEnd` works under Wolfram Engine 14.3 in this container — that is **one live
check**, and if it fails the fix is to guard the PDF export rather than to abandon the run:
the analysis and the association are what Stage 1 needs, and the PDF is a human
cross-check.

---

## 1. Environment audit (verified 2026-09-03)

| component | state | note |
| --- | --- | --- |
| Wolfram Engine | **14.3.0 installed and activated** | `wolframscript` 1.13.0; licensing present under both `~/.Mathematica` and `~/.WolframEngine`. PSALTer requires 14.0+ |
| xAct | **1.2.1 installed** in the Wolfram user base | PSALTer requires 1.2.0+; xTensor, xPerm, xCoba, xPert, xTras, SymManipulator all present |
| PSALTer | **not installed** | v2.0.2, GPL-3.0-or-later |
| Inkscape | **not installed** | not actually needed (§0.6) |
| `sudo -n` | works | apt install is feasible |
| network | GitHub reachable | full clone of `SupplementalMaterials-2607` is slow (large PDF/blob payload) — prefer raw-file fetches or sparse checkout |
| `.wxf` fixtures | **not in repo** | upstream only; 692 B and 2,874 B, sizes match H6's decode notes |
| H6 decode scripts | **salvaged** | were in `/tmp`; now `scripts/research/psalter_stage1/wxf_decode/` |

**Third-party discipline.** PSALTer and the supplemental materials are GPL-3.0-or-later and
TIDAL is MIT; the trigger for the license question is *distribution*. We therefore commit
**the route, not the payload**: `scripts/research/psalter_stage1/fetch_reference_sources.sh`
pins the revisions and downloads into a gitignored working directory. Reading and adapting
with provenance is authorized (D6); redistributing is the release blocker recorded on #495.

## 2. Install route

**`scripts/install-psalter.sh`** — model it on `scripts/install-xact-xcoba.sh`, which
already has the shape: colored log helpers, a `check_wolfram` guard, and a
`get_wolfram_user_dir` that *queries* `$UserBaseDirectory` through `wolframscript` rather
than hardcoding a path (required — no environment-specific absolute paths in committed
files).

1. `git clone` PSALTer at a **pinned commit** (default `bb45adb0`, overridable by
   environment variable), into a temporary directory.
2. Copy `xAct/PSALTer` into `$UserBaseDirectory/Applications/xAct/` — that is the entire
   install, per the upstream README (its own instructions show the global variant,
   `sudo cp -r PSALTer /usr/share/Wolfram/Applications/xAct/`; the user base is the
   no-sudo equivalent and matches how xAct is already installed here).
3. Write the resolved commit hash to `INSTALLED_COMMIT` inside the installed directory —
   the exporter reaches into private symbols (§0.4), so knowing exactly which revision is
   installed is a correctness input, not bookkeeping.
4. `sudo apt-get install -y --no-install-recommends inkscape`, warn-and-continue.
5. Verify: `Needs["xAct`PSALTer`"]` loads, and `DefField` and `ParticleSpectrum` both
   exist. Allow a generous timeout — the first load builds SPO tables and is slow.

**`scripts/verify-wolfram-setup.sh`** gains a PSALTer section doing step 5 plus a
`DefField` smoke test on a scalar, so the check survives independently of the installer.

## 3. Install verification — three tiers, cheapest trust first

The point of this task is that a *partially* working install must fail loudly. H6 chose
TorC because one of its three conditions is a **massless residue** condition, so an install
whose massless analysis is broken cannot pass quietly. Preserve that property in every
tier.

**Tier 1 is the install gate** (settled 2026-09-04). The reason is what makes the choice
non-arbitrary rather than merely cheapest-first: *an install gate must fail unambiguously.*
Tier 1's input and its expected result are both the author's, so a mismatch can only be the
install. Tier 2 mixes install correctness with our own authoring, so a Tier-2 failure does
not localize — the same "never conflate two hypotheses in one test" principle that makes the
frozen legacy specs usable as physics oracles. **Tier 2 is therefore the physics gate, not
the install gate**, and Tier 3 substitutes for Tier 2 when TorC's ECT authoring is deferred.

**Tier 1 — replay the published input (cheapest, strongest).** Run
`ParticleSpectrographCTEG.m` unmodified (fetched, not committed), then load the
repository's committed `ParticleSpectrographCTEG.mx` and diff the two associations
key-by-key. This tests the install against the author's own artifact with zero physics
authored by us. It is the tetrad/PGT formulation (21 generators).

**Tier 2 — the ECT reproduction (our formulation, the design's primary oracle).** Author
the post-Riemannian run — fields `h`, `T`, `a`; 15 gauge generators — using the §0.5
single-session pattern, modeled on 2506b's `PoincareGaugeTheory.m` + `Linearise.m`, with
the coupling dictionary that makes `eq:TorCLagrangian` linear per operator (Barker's own
`MuLambda` is the precedent). Four machine-checkable assertions:

1. **gauge generators = 15** — the count of source-constraint rows
   (`ConstructSourceConstraints` → `SourceConstraintRows`);
2. **massive content = one `0⁻` state**, with `m²` mapping to `λ` under the dictionary
   (`Simplify` to exact equality, not numeric agreement);
3. **massless polarizations = 2**;
4. **unitarity conditions logically equivalent** to `λ ≥ 0 ∧ μ < 0 ∧ (ν+2μ)(ν−μ) > 0`
   after dictionary substitution, via `Reduce[exported ⇔ published, vars, Reals]` under
   `MPlanck2 > 0`. Robust to algebraic form, and it is the assertion that carries the
   massless condition — the one that makes a half-working install fail loudly.

**Tier 3 — the fully local alternative.** `literature/2506.02111/` carries
`ScalarParityViolatingPGT.tex` (a literal `ParticleSpectrum[…]` listing),
`TetradPerturbation.tex` + `SpinConnection.tex` (the `DefField` declarations),
`LoadPSALTer.tex`, and `ParticleSpectrographScalarParityViolatingPGT.pdf`. Its published
massless condition is `c₁ < 0`, so it keeps the fail-loudly property with zero fetching,
and H6's addendum records it as an **equal option, not a degraded substitute**. One trap:
its `SpinConnection` is declared `Antisymmetric[{-b,-c}]` while CTEG's is
`Antisymmetric[{-a,-b}]` — copy each listing exactly as published rather than normalizing
across them.

**The PDF is a recorded eyeball check, not a gate.** Comparing our rendered spectrograph
against `literature/2507.09228/ParticleSpectrographTorCECT.pdf` is worth doing and worth
recording; programmatic PDF comparison is not realistic and the assertions above are what
must pass.

## 4. The spectrum derivation branch

### 4.1 Where it lives

`tidalcosmo/` (README stubs today; never build on the legacy layout). Specifically:

- **`tidalcosmo/config/`** — the *shared, convention-free* input model of H6 §4.5: field
  declarations and symmetry classification, coupling roster, term-level Lagrangian
  structure, scope guards. Frozen dataclasses, per H4 §2.8's PSALTer-derived idiom. Mark
  every docstring **provisional pending WS1/M0** — this is a minimal stand-in for a typed
  config layer that WS1 owns, deliberately built now so H8 does not block on M0.

  **Build it, but mark it disposable** (settled 2026-09-04). WS1 is next in dependency order
  and not dispatched, so "imminent" is not knowable and waiting could stall indefinitely. The
  risk worth naming out loud is the one that actually happens to provisional layers: they
  become the de facto config **by accretion**, one convenience field at a time, until
  replacing them is a migration. So the requirement is explicit and testable — **WS1's config
  surface replaces this module; it does not extend it.** Keep it minimal enough that
  replacement is deletion: no field that Stage 1 does not need, no downstream module importing
  it except the spectrum branch, and the marker in every docstring so a later reader cannot
  mistake it for the design.
- **`tidalcosmo/derive/`** — the generator, the `wolframscript` driver, and the committed
  `.wl` exporter package.
- **`tidalcosmo/spectrum/`** — the Python side of the contract: the WXF reader and the
  contract dataclasses that Stage 2 will consume.

`pyproject.toml`'s `include = ["tidal*"]` already matches `tidalcosmo`; the only packaging
change is adding the `.wl` files to `package-data`.

### 4.2 The input model, and the reject rule

Term-level, not one expression string: `[[spectrum.fields]]`, `[spectrum.couplings]`,
`[[spectrum.lagrangian.terms]]` with each term naming its coupling and its operator.
Validate at load with `error_with_hint`:

- **Field symmetry classes are an enumeration, and it is now known exactly.** PSALTer
  registers rank 0; rank 1; rank 2 general/antisymmetric/symmetric; rank 3 general,
  `Antisymmetric{12,13,23}`, `Symmetric{12,13,23}`, totally antisymmetric, totally
  symmetric (`Sources/DefField.m:11-24`). Torsion is `RegisterFieldRank3Antisymmetric23`.
  Validate against exactly this list and name the legal values in the hint.
- **Coupling per operator; reject, never auto-assign** (H6 §4.1). Since PSALTer will not do
  this for us (§0.2), our validator is the only guard. A term with a bare numeric
  coefficient is a config error naming the term, with a hint explaining that PSALTer's
  theory space is `S = ∫ Σᵢ θᵢ Oᵢ` and that the numeric value belongs in the coupling's
  *value*, not the Lagrangian.
- The legacy `examples/torsion_gertsenshtein/theory.toml` is a **negative-test fixture**:
  its `(1/kappa^2) R̃ … − ¼F²` string must be rejected. The roster-conformant restatement
  is a new `theory_spectrum.toml` — Einstein–Hilbert with its own generic coupling (H6
  §4.1: no distinguished κ), `α₁I1 + α₂I2 + α₃I3`, and `cF·F²`, **without b5**.

### 4.3 The emitted script

Single session (§0.5), standalone, in Barker-native conventions from its first line —
signature `(+,−,−,−)`, `ε₀₁₂₃ = +1` — which cost nothing because PSALTer establishes that
geometry itself. Emit, in order:

1. `Needs["xAct`PSALTer`"]`, then `DefConstantSymbol` per declared coupling and `DefField`
   per declared field (with `PrintAs`/`PrintSourceAs`, so the spectrograph is readable).
2. The nonlinear geometric objects as `DefTensor`s on PSALTer's `M4`, with `MakeRule`
   definitions in terms of the perturbation fields — the ported analogue of
   `tidal/cli/_derive.py::_wls_torsion_curvature_decomposition` (line ~2088,
   `ChangeCurvature[L, CDT, CD]` plus the contortion identity
   `K^a{}_{bc} = ½(T^a{}_{bc} + T_b{}^a{}_c − T_{bc}{}^a)`). **Ported, never imported**
   (D3), with a provenance docstring citing the legacy function; its physics comments carry
   over, its `(−,+,+,+)` sign conventions do not.
3. Bare `h_{μν} ≡ g_{μν} − η_{μν}`, plus sign, no rescaling (H6 §4.2); derived fields
   pre-expanded (`F = dA`).
4. Linearization: multiply the perturbation fields by a bookkeeping constant,
   `Series[…, {ε, 0, 2}] // Normal`, set `ε → 1`, canonicalize.
5. **Scope guards, in-script**: assert nothing linear in the perturbation fields survives;
   assert no background-field symbol appears (the spectrum branch consumes the vacuum
   Lagrangian block only); assert the declared signature and `ε` orientation, printing both
   into the artifact header so the convention is a machine-checked property of the run and
   not a comment.
6. `ParticleSpectrum[…, TheoryName→…, MaxLaurentDepth→1, Method→"Hard"]` wrapped in
   `AbsoluteTiming` with timestamped checkpoints, then `Get` our exporter.

The gauge-unfixed Lagrangian is handed over deliberately: finding the gauge symmetries and
imposing source constraints is PSALTer's job (H6 §4.5).

### 4.4 The driver

`tidalcosmo/derive/wolfram_driver.py`: engine-idle guard before every launch
(`pgrep -f -i 'wolframscript|WolframKernel|MathKernel'` — the single-license rule is
project-wide, not per-session), per-run log capture, wall-clock timing, configurable
timeout, strictly serial. Ported in shape from `tidal/cli/_derive.py::_run_wolframscript`
(line ~6438). No parallelism anywhere in this subsystem, ever.

## 5. The exporter

**`tidalcosmo/derive/wolfram/Stage1Export.wl`**, with a provenance header recording that
the structure-harvesting pattern is modeled on `SupplementalMaterials-2607`
`WolframLanguage/ParticleSpectroscopy/JuliaExport.m`, and that the coefficient-tensor and
explicit-label export exist nowhere upstream.

**Input surface** (§0.4): the theory association (`<Name>@WaveOperator`,
`@PseudoDeterminant`) plus the `xAct`PSALTer`Private`` globals `$LocalSourceConstraints`,
`$LocalMasslessSpectrum`, `$LocalSpectrum`, `$LocalUnresolvedPoles`,
`$LocalOverallUnitarity`, `$LocalSummaryOfTheory`. Read defensively — normalize keys, treat
a bare unevaluated symbol as *absent*, and handle both the populated case and the
degenerate one (`UnitarityConditions = Text["(Demonstrably impossible)"]` with an empty
constraint matrix is a legal output).

**Output: WXF**, one file per theory. Exact rationals survive, `wolframclient` decodes it,
and it is the format Stage 2's ingest already speaks. Contents, per H6 §6.1:

- `schema_version`, `theory`, `psalter_commit`, `conventions`
  (`signature`, `epsilon0123`), ordered `couplings`, `wall_time_s`, per-stage checkpoints.
- `sectors`: for each, the states with **explicit per-state `J^P` labels** (never
  positional), the coupling-linear coefficient tensor
  `A[(n_couplings+1) × 3 × dim × dim]` such that the block is
  `A[0] + Σᵢ cᵢ·A[i]` with the second axis indexing powers of `k`, the gauge rank, and the
  massive/massless index partition that the Schur reordering needs.
- `source_constraints` (the constraint matrix; verify coupling-independence and flag if
  not — H6 §9 relies on structural gauge symmetries being coupling-independent),
  `gauge_generators`, `massless_polarizations`, `massive_spectrum`,
  `unitarity_conditions` as a status enum plus expression, per-sector `z_degree`
  (`Exponent` in `k` — cheap symbolically at Stage 1; document the non-generic-point
  failure mode), and `k4_guard`.
> **Amendment (coherence pass, 2026-09-04).** The `A[0]`-vanishes acceptance gate below is
> conditional, and the condition was never stated. If every term carries exactly one
> *sampled* coupling then the coupling-free slot must be zero, and a non-zero `A[0]` is
> precisely the bare `−¼F²` that PSALTer passes silently (#522) — the gate is then the
> replacement for the check PSALTer does not perform. **But TorC-class Lagrangians carry
> terms in `M_p` and `Λ`**, which are theory-fixed constants rather than sampled couplings;
> if those sit outside the coupling vector they land in `A[0]` *legitimately* and the gate
> is wrong as written.
>
> **Recommendation:** non-sampled constants stay OUT of the coupling vector (they are fixed,
> not sampled), so `A[0]` is legitimate and the coupling-linearity check becomes *"no sampled
> coupling appears non-linearly"* — a different test. **Confirm against a real export before
> relying on either form**; this is the one item the documents alone cannot settle. Matching
> amendment in `spectrum_design.md` §6.1.

- **An in-kernel reconstruction assertion**: reassemble `A[0] + Σᵢ cᵢ·A[i]` and `SameQ`
  it against the original block after `Together`. An exporter that silently drops a term is
  the failure mode that would poison every downstream verdict.

**`J^P` labels — the plan for H6's open item §14.2.** The label is not decoration: parity
feeds the Hermitization `diag(1,−1)` and hence every ghost verdict, and the released JAX
`SPIN_LABELS = {"0+","1-","2+"}` positional convention is *demonstrably* wrong on the
release's own `A23Theory` fixture. `JuliaExport.m` confirms the defect verbatim in its
comments ("Index 1 = spin-0+"). Four steps, in order:

1. Dump PSALTer's field-kinematics structures for `h`, `T`, `a` on the live install and
   compare against the published `FieldKinematics*.pdf` tables in `literature/2507.09228/`
   and `literature/2506.02111/`.
2. Read the labels from the SPO/field-kinematics tables (private symbols; safe because the
   commit is pinned, and recorded as a dependency).
3. **Calibration tests** — the step that turns the convention from assumed to measured.
   Run single-operator theories whose block occupancy is known a priori: a pure mass term
   `Θ₃A_μA^μ` (occupies `0⁺` and `1⁻` only), a pure torsion trace term, and the design §8
   Vector theory whose two 1×1 blocks are published verbatim
   (`½[k²(Θ₂−Θ₁) − Θ₃]` and `½[−k²Θ₁ − Θ₃]`). The labeled rows must be exactly the
   non-zero ones.
4. Reproduce the `A23Theory` fixture's three `J`-blocks of dims 2/4/2
   (`2·1 + 4·3 + 2·5 = 24`) and label the `J=1` block's visible `2+2` split as `1⁺ ⊕ 1⁻`.

**Fixtures.** Commit the two small upstream `.wxf` files under `tests_cosmo/fixtures/psalter/`
with a `PROVENANCE.md` (source repository, pinned revision, file hashes, D6 note, and the
observation that these are *data outputs*, not GPL sources — flagged on #495 for the
release-time review regardless). They are answer keys that let the Python reader be tested
in CI **without Wolfram installed**, which is the property that matters.

## 6. Cost measurement protocol

Read §0.1 first: the `Method` framing is wrong on the pinned release, and the measurement
must be reported for what it is.

- **Subject:** the roster-conformant PGT+EM theory — `theory_spectrum.toml`, EH + `α₁I1 +
  α₂I2 + α₃I3` + `cF·F²`, **no b5** — generated end to end by our own branch.
- **Instrumentation:** `AbsoluteTiming` around the whole call plus timestamped checkpoints
  at PSALTer's own stage boundaries (`ConstructWaveOperator`, `ConstructSourceConstraints`,
  `ConstructSaturatedPropagator`, `ConstructMassiveAnalysis`, `ConstructMasslessAnalysis`,
  `ConstructUnitarityConditions`, `ConstructSpectrograph`). Partial timings localize the
  cost even if the run never finishes, which is why they matter more than the total.
- **Scheduling:** background, engine-idle guard first, nothing else touching Wolfram.
- **Ceiling:** review at ~1 h and ~8 h; hard abort at **24 h**. If it does not terminate,
  report "did not terminate within 24 h" plainly with the checkpoint trace. Do not work
  around it, do not silently shrink the theory — the design owns that contingency and the
  honest number is the deliverable.
- **Free extra data point:** the Tier-1 and Tier-2 oracle runs (§3) are PGT-sized and
  include curvature-squared terms; time them too and report all of them together.

**Missing-coupling probe** (minutes, after the main run): feed PSALTer a term with a bare
numeric coefficient, bypassing our validator. Expected outcome per §0.2 is **silence** —
record that, and record the `UnknownCoupling` text for a genuinely undeclared symbol, which
is the message our own hint should echo.

**The b5 probe is a guard test, not roster support.** Curvature-squared terms are excluded
from the roster by user decision (`spectrum_design.md` §2, lines 83-87). Run
`theory_spectrum_b5.toml` once, bounded at 8 h, and record the outcome as **"guard fires
correctly"**: either PSALTer produces `k⁴` entries and our `k4_guard` flags the sector, or
PSALTer errors and we transcribe it. Both are passes. **Do not conclude that b5 is
supported because the test exercised it** — that inversion is the specific failure this
addendum exists to prevent.

**Where the numbers land:** a new `docs/cosmology/stage1_measurements.md` — host CPU/RAM,
Wolfram and xAct versions, PSALTer commit, every wall time with its checkpoint trace, the
failure-mode transcripts, and one explicit go/no-go sentence against the design's standard
("minutes are fine, hours are survivable"). `CAMPAIGN.md` is for campaigns and is
append-only; a handoff report is ephemeral; this needs to be a committed document.

## 7. Gates, stated before any code is written

| gate | criterion |
| --- | --- |
| install | committed script produces a clean `Needs["xAct`PSALTer`"]` on a fresh user base; version banner printed |
| oracle T1 | our `CTEG` run's association matches the committed `.mx` key-by-key |
| oracle T2/T3 | gauge-generator count exact; massive content with dictionary-exact `m²`; massless polarization count; `Reduce`-verified equivalence of unitarity conditions **including the massless one** |
| reader | `A23` blocks are dims (2,4,2) with `2·1+4·3+2·5 = 24`; `Vector` blocks equal the published expressions exactly (symbolic difference simplifies to zero); placeholder/plural-key/degenerate cases unit-tested |
| reject rule | a bare-numeric term errors with a hint naming the term; the legacy `theory.toml` Lagrangian is rejected |
| exporter | in-kernel reconstruction `SameQ`; Python-side numeric agreement ≤ 1e-12 at 5 random rational coupling points; label calibration passes on ≥ 3 single-operator probes; `A[0]` vanishes **iff non-sampled constants are excluded from the coupling vector** — see the amendment below (coupling-linearity) |
| emitted script | declares and machine-checks signature and `ε`; standalone; no repo-absolute paths |
| cost | measured, not guessed, with checkpoints; one explicit go/no-go sentence; `Method` inertness confirmed live |

**Where tests live.** `scripts/run_wolfram_tests.sh` auto-discovers `tests/wolfram/test_*.wls`
by glob, so a new suite is picked up automatically — keep it **fast** (Vector-scale only,
seconds to minutes) or it will slow every future full-suite run. The oracle and PGT+EM runs
belong in `scripts/psalter/` wrappers, outside the glob. Python tests skip cleanly when
`wolframclient` is absent.

## 8. Risks

| risk | mitigation |
| --- | --- |
| `UsingFrontEnd@Export` fails headless | one live check; if it fails, guard the PDF export — the association is what Stage 1 needs (§0.6) |
| private-symbol churn in PSALTer | pin the commit, record it in every export, and let the fixture tests fail loudly on drift (§0.4) |
| single-session linearization proves awkward for our field content | documented two-session fallback (§0.5); last resort, hand-derive the (small) quadratic Lagrangian once and automate only the spectroscopy — report the compromise rather than hiding it |
| `ParticleSpectrum` does not terminate on PGT+EM | report plainly with checkpoints; no workaround (§6) |
| `Method` turns out to be live in a newer revision | §0.1 step 2 confirms on the install; re-time if so |
| license | PSALTer never committed; fetch script commits the route; `.wxf` data fixtures and adapted `.wl` code carry provenance; #495 comment for the release review |
| one-license contention | engine-idle guard in the driver; long runs backgrounded and serialized |

## 9. How this reaches the sampler

Stage 1 exists to make Stage 2 cheap, and Stage 2 exists to be **a pre-run validity check
inside the sampling loop**: a coupling point that is vacuum-sick is rejected before any
forward solve is spent on it. That is H6 §11 ("screen before compute") and §2's one-way
filter — vacuum-sick is a reason to discard; vacuum-healthy is not a certificate about
behavior on an FRW background. Mechanically the verdict is a flagged rejection through the
same mechanism `validity/` already uses, carrying its reason, and it composes with the
soft-likelihood option (H6 decision 6) if the program later wants unitarity to shape the
posterior rather than gate it. Everything in this document is in service of that check
being correct and taking a millisecond.

## 10. Recommended step sequence for the implementing session

The critical path is Wolfram-serial; Python work fills the gaps.

1. **Install** (§2) → commit `feat(scripts): add PSALTer install and reference-fetch scripts`.
2. **Tier-1 oracle** (§3), backgrounded. While it runs: **fixtures + WXF reader + contract
   dataclasses + pytest** (§5) → commits `test(spectrum): …` and `feat(spectrum): …`.
3. **Tier-2 oracle** authored and run. While it runs: **input model + validator + generator
   + driver**, text-emission only, with golden-text tests (§4) → commit
   `feat(config,derive): …`.
4. **Exporter** developed against Vector-scale theories, serial short runs (§5) → commit
   `feat(derive): …` with the Wolfram test suite entry.
5. **PGT+EM cost run**, backgrounded (§6). While it runs: docs, provenance,
   measurements-document skeleton.
6. **Missing-coupling probe**, then the **b5 guard probe** (§6) — both after the main run,
   both short.
7. **`docs(cosmology): record Stage-1 measurements and the cost verdict`**, plus README
   updates for `tidalcosmo/derive/` and `tidalcosmo/spectrum/`.

The three findings that were actionable at study time are already filed — #521, #522 and
#523 — so the implementing session inherits them rather than rediscovering them. File
further issues as they appear, and add the `.wxf`-fixture license note to #495 when the
fixtures are committed. Report the branch to the orchestrator; do not merge.

## 11. Decisions (settled by the user, 2026-09-04)

The study closed with three open questions. All three are answered; the reasoning is recorded
because in each case it generalizes past this handoff.

1. **Amend the design, pinned to the version.** The `Method` finding lives in
   `spectrum_design.md` §4.5/§14.1 — not only in the measurements record. *An observation can
   live in a record; a finding that invalidates a stated instruction has to live where the
   instruction is*, because §4.5 is what an implementer reads. Pinned to `bb45adb0` because a
   later PSALTer release could implement the option, and a bare "it's dead" would then be
   wrong in the other direction. Applied 2026-09-04: six inline amendment notes covering
   §4.1 (#522), §4.5 (#521), §6.1 (#523) and §14.1/§14.2. The same principle promoted the
   other two findings from "recorded" to "amended" — each of them also invalidated an
   instruction, not merely a background claim.
2. **Tier 1 is the install gate** (§3). An install gate must fail unambiguously; Tier 1's
   input and result are both the author's, so a mismatch can only be the install, while a
   Tier-2 failure is ambiguous between install and authoring. Tier 2 is the physics gate.
3. **Build the provisional config layer, and mark it disposable** (§4.1). WS1 is not
   dispatched, so blocking on it could stall indefinitely; the named risk is a provisional
   layer becoming the de facto config by accretion, and the guard is that WS1's surface
   **replaces** it rather than extending it.

Two program-level decisions remain with the user and are unaffected by this study: the rung
order after O1, and when to start the `tidalcosmo` version line.
