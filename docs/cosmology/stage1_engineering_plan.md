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

`ValidateLagrangian.m` defines six messages and has throw sites for four — but only
**three can fire on this install**, per the amendment below:

| message | fires on this install? | throw site in the source? | fires on |
| --- | --- | --- | --- |
| `Zero` | yes | yes | `PossibleZeroQ` on the Lagrangian |
| `UnknownCoupling` | yes | yes | a symbol that is not a `ConstantSymbolQ` and not an `xTensorQ` |
| `UnknownField` | yes | yes | a tensor with no `xAct`PSALTer`<name>`` context — i.e. not `DefField`ed |
| `NonQuadraticFields` | **no** — guard unavailable, see below | yes | `PolynomialDegree > 2` in the fields |
| `NonLinearCouplings` | no | **never** | — (message defined at line 6, no throw site in the package) |
| `ParityOdd` | no | **never** | — |

The first column is the one that matters when reading this table on a live install, so it
comes first; the second records what the *package* would do given its dependencies.

> **⚠ Amendment (I-REM, 2026-09-09 — PSALTer v2.0.2 `bb45adb0`).
> `EXPIRES-WITH: #551` → EXPIRED 2026-09-11: in the certified configuration this guard
> FIRES.** `scripts/psalter/register_resources.wl` supplies the resource locally, and a cubic
> control now throws where an unregistered install accepted it. What follows describes an
> **unregistered** install — still the state of any fresh container until the registration is
> re-run, which is why `install-psalter.sh` now does it. `NonQuadraticFields` has a throw site,
> but it is guarded by
> `ResourceFunction["PolynomialDegree"]` (`ValidateLagrangian.m:38`), and the Wolfram
> Function Repository cannot be reached from this container (#551). PSALTer emits
> `ResourceObject::notfname` and **carries on**, so the check silently never runs. The
> distinction that matters for anyone reading this table: *the source has four throw sites;
> this environment fires three.* The fourth returns the moment the resource becomes
> available. Counted with the other silently-absent validations in
> `stage1_measurements.md` §6.

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

> **⚠ Amendment (I-REM, 2026-09-09) — this is the one §0.x finding Wave 0 did not test, and
> it now has an owner.** `stage1_measurements.md` §8's scorecard records every §0.1–§0.6
> finding as confirmed or overturned against the live install, except this one: "both
> linearization architectures are documented — **not tested**, out of scope for Wave 0". It
> remains a static source reading. **Owner: I-S1B** (the Stage-1 Wolfram-side handoff,
> #495), which is the first session that will actually exercise a linearization
> architecture and can therefore settle it as a by-product rather than as separate work.

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
  `xAct`PSALTer`Private`ToNewCanonical`. The README's own words (quoted verbatim, British
  spelling and all — the project writes American English, but quoted sources are not
  rewritten):
  <!-- cspell:disable -->
  "one can always obtain the correct linearisation by reinterpreting the model as a
  field theory on flat spacetime:
  this approach is guaranteed to work because particle physics is not actually sensitive to
  geometry per se."
- **Two sessions, copy the expression across**: also documented — "prepare the linearised
  expression in a separate *xAct* session, taking advantage of the full geometric
  interpretation, and then copy the resulting *Wolfram Language* expression directly into
  the PSALTer session. Care must be taken in this case to ensure that the indices are
  correctly matched."
  <!-- cspell:enable -->

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

> **⚠ Amendment (I-526 live check, 2026-09-06 — PSALTer `bb45adb0`, Wolfram 14.3, this
> container).** The live check was run and **the prediction above is wrong in the way that
> matters.** `UsingFrontEnd@Export` does not *fail* headless — it **hangs indefinitely**,
> and the fix is **one environment variable**, not a guard:
>
> ```sh
> export QT_QPA_PLATFORM=offscreen   # export then completes in ~4 s
> ```
>
> ~~`DISPLAY` **is** set in this container (`:27`), so Qt finds a display, tries to use it,
> and blocks.~~ **That mechanism is wrong — see (b) below.** It was written here by the
> orchestrator from correlation (the variable is set) without being tested, one day before
> the programme adopted the rule that anything written into an instruction site is an
> assertion and gets verified first. The fix and its effect were right; the cause was not.
>
> **Why the distinction matters more than it looks.** "Guard the PDF export" would have had
> the next session *delete a working capability* to avoid a problem that has a one-line
> fix. And a hang, unlike a failure, consumes the **single-license Wolfram lane
> indefinitely** — inside the §6 cost run, with its 24 h ceiling, an unguarded hang is
> indistinguishable from "PSALTer is slow on PGT+EM", which is precisely the measurement
> that run exists to make.
>
> **Extended (I-526 implementation, 2026-09-07).** Two corrections to the
> amendment above, both from the live install.
>
> **(a) The exposure starts at `DefField`, not at the spectrograph.** `DefField`
> calls `SummariseField` (`Sources/DefField.m:91`), which exports a
> `FieldKinematics<Field>.pdf` through the front end
> (`Sources/DefField/SummariseField.m:85`), unguarded. So `ParticleSpectrographCTEG.m`
> would hang at **line 7 of 10**, seconds in — not hours later. Good news, in that
> it fails fast; but the variable is needed far more broadly than "before running
> `ParticleSpectrum`".
>
> **(b) `DISPLAY` is not the cause.** The cause is that no Qt *platform plugin*
> could be initialized: `offscreen` was the only one with all its shared libraries
> present (`xcb` is missing five, `wayland-egl` was missing `libwayland-egl1`).
> Proved by forcing `QT_QPA_PLATFORM=xcb`, which reproduces the abort-then-block
> exactly. With the libraries present the same `DISPLAY=:27` works in ~1.2 s.
> Note the optional Inkscape step installs `libwayland-egl1` and so fixes this by
> side effect — luck, not a guarantee; keep setting the variable.
>
> Measurements: `stage1_measurements.md` §2.3–2.4.

> **Where the export lives:** in the scripts (`install-psalter.sh`, the Tier-1 gate script,
> any PSALTer runner), *not* in `devcontainer.json` — a container-wide setting would
> override a display other processes may use, and would need a rebuild to take effect, so
> it would look applied while being inert. Durable coverage belongs in
> `tidalcosmo/derive/wolfram_driver.py` (§4.4), beside the engine-idle guard.

---

## 1. Environment audit (verified 2026-09-03)

| component | state | note |
| --- | --- | --- |
| Wolfram Engine | **14.3.0 installed and activated** | `wolframscript` 1.13.0; licensing present under both `~/.Mathematica` and `~/.WolframEngine`. PSALTer requires 14.0+ |
| xAct | **1.3.0 bundle installed** in the Wolfram user base (corrected 2026-09-11 — the "1.2.1" audited here was the install script's default, not a measurement; xTensor 1.3.0 / xPerm 1.2.4 / xCore 0.6.10 / xCoba 0.8.6) | PSALTer requires 1.2.0+; xTensor, xPerm, xCoba, xPert, xTras, SymManipulator all present |
| PSALTer | **not installed** | v2.0.2, GPL-3.0-or-later |
| Inkscape | **not installed** | not actually needed (§0.6) |
| `sudo -n` | works | apt install is feasible |
| network | GitHub reachable | full clone of `SupplementalMaterials-2607` is slow (large PDF/blob payload) — prefer raw-file fetches or sparse checkout |
| `.wxf` fixtures | **not in repo** | upstream only; 692 B and 2,874 B, sizes match H6's decode notes |
| H6 decode scripts | **salvaged** | were in `/tmp`; now `scripts/research/psalter_stage1/wxf_decode/` |

> **⚠ Amendment (I-REM, 2026-09-09) — three rows above are superseded, and two facts are
> added. The table itself is left as written: it is a dated audit ("verified 2026-09-03")
> and rewriting its rows would destroy the record of what was true when this plan was
> drafted.**
>
> **Superseded rows**, all re-verified 2026-09-09:
>
> | row | as audited | now |
> | --- | --- | --- |
> | PSALTer | not installed | **v2.0.2 `bb45adb0`, installed and verified 2026-09-07** (#526) — **certified 2026-09-11** (I-543, #557): Tier 1 `VERDICT: MATCH` once the two undeclared Function Repository resources are registered locally; the mismatch this row recorded on 09-09 is the failing evidence at `evidence/tier1-20260907/` |
> | Inkscape | not installed | **1.2.2-2+b1**, at `/usr/bin/inkscape` |
> | `.wxf` fixtures | not in repo | **two are committed**, at `tests_cosmo/fixtures/psalter/` (`A23Theory`, `VectorTheory`), with `PROVENANCE.md`; flagged for the release-time license review on #495 |
>
> **Added — the Wolfram Cloud is unreachable from here** (`EXPIRES-WITH: #551` → **still
> true 2026-09-11, and NOT what was fixed**). The site itself recovered — it answered 503 on
> 09-09 and 200/302 on 09-11 — but the repository served **no definitions to anyone** during
> I-543's session, including to an authenticated login and to in-cloud evaluation. #543 was
> resolved by registering the genuine definitions locally, **not** by restoring reachability,
> so anything needing a *new* resource still hits this. Measured 09-09: 503 from
> `www.wolframcloud.com` and 404 from the resource API, from **both `curl` and Wolfram**,
> while GitHub returns 200 — so this is network-layer, **not** Wolfram authentication.
> Consequence: `ResourceFunction[…]` downloads fail, PSALTer emits `ResourceObject::notfname`
> and **carries on**, and the guarded checks silently never run (§0.2's table, and the
> `LinearlyIndependent` sites under `ConstructSourceConstraints/`). **Function Repository
> resources must therefore be supplied locally**. I-526 reported doing so on a throwaway copy
> and getting an identical verdict, but recorded only *that* it supplied them, never *how* —
> and the sites are reached partly through `NewParallelSubmit`, so a definition present only
> on the master kernel would not exist on the subkernels. **Whether that stand-in was ever
> live was re-verified by I-543 and the refutation was **invalid**: the stand-in reproduced the disabled behavior (#556, closed), so the missing resources were the cause all along (#543, closed 2026-09-11).
> Why it is unreachable is **not diagnosed** — it is not established that it is impossible, only
> that it does not work as configured — so this is tracked as **#551** rather than treated as
> a standing constraint.
>
> **Added — Inkscape's *dependencies* are what make the headless export work**, which the
> "not actually needed (§0.6)" note does not convey. Installing it pulls in
> **`libwayland-egl1`**, which makes the `wayland-egl` Qt platform plugin loadable and fixes
> the export hang by accident. §0.6 is right that PSALTer's own `$InkscapePath` is unused,
> but `--skip-inkscape`, a slimmer image or another distribution reintroduces the hang.
> `QT_QPA_PLATFORM=offscreen` (§4.4) is the fix to rely on; Inkscape is not.

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

> **⚠ Amendment (I-REM, 2026-09-09) — what a Tier-1 pass would certify, and what it
> would not. This half carries NO expiry: it is a permanent property of a two-key
> oracle, true whether the gate passes or fails.**
>
> Call it the **Tier-1 `.mx` oracle** rather than "the oracle" — this document and
> `repo_reshape.md` §8 use that word for two unrelated artifacts, the other being the
> **frozen legacy oracle** of the M0.5 milestone (185 committed files of `tidal` output).
>
> `ParticleSpectrographCTEG.mx` holds exactly **two** keys — `WaveOperator` and
> `PseudoDeterminant` (§0.4, measured). Everything else `ParticleSpectrum` computes lives in
> private globals the `.mx` never contains. So **a passing Tier 1 certifies field
> declaration, decomposition, wave-operator construction and the pseudo-determinants, and
> nothing else** — not source constraints, not the spectrum, not the unitarity conditions.
> Tier 1 was never going to certify as much as its name suggests, which is why §3 calls
> **Tier 2 the physics gate**.
>
> **Consequence, and a requirement on Tier 2:** the 0-versus-21 source-constraint defect
> (#543) sits in exactly that uncovered region, and no fixture would have caught it — it was
> found by reading the private global directly and comparing against a count in the
> *published formulation* (the papers and PSALTer's own README and supplemental material),
> which is human-readable published source rather than a byte-comparable artifact.
> **Tier 2 must therefore add source-constraint coverage.** Its assertion 1 below already
> names the mechanism — gauge generators = the count of `ConstructSourceConstraints` →
> `SourceConstraintRows` — so this is a requirement to keep it, not to invent it.
>
> **Separately, and this half DOES expire — `EXPIRES-WITH: #543` → EXPIRED 2026-09-11:** Tier 1
> ~~currently **fails** on this install~~ **now passes** (`MATCH`, both keys identical) on
> Wolfram 14.3.0 × PSALTer v2.0.2 `bb45adb0`, with PSALTer's two undocumented Function
> Repository dependencies registered from genuine code (`scripts/psalter/register_resources.wl`).
> The 2026-09-07/09 failure — `WaveOperator` bit-exact, `PseudoDeterminant` all zeros, reproduced
> independently — was those missing dependencies, not the engine; the gate was never relaxed and
> no workaround was applied. Record: `stage1_measurements.md` §4.4/§4.8, evidence
> `docs/cosmology/evidence/tier1-20260911-pass/`, `psalter_543_investigation.md`.

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

> **⚠ Amendment (orchestrator, 2026-09-15 — D-A, adopted from R-1 #566): nothing is emitted.**
> The six steps below remain the order of work, but they are performed by **one committed
> Wolfram package** (`RunStage1`) run by a **fixed `driver.wls`**, receiving the theory as WXF
> data. Names, the symbol whitelist and coupling-linearity are checked **in the package, on a
> held parse, before PSALTer runs** (`interfaces_decision.md` §2.4 — eleven malformed and
> injected inputs refused). Step 5's signature/`ε` assertion becomes the spectrum manifest's
> `conventions` record, asserted by every consumer (`conventions.md` §5). Read "emit" below as
> "the package performs".

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

> **⚠ Amendment (orchestrator, 2026-09-15 — D-A).** The driver launches the fixed
> `driver.wls <theory.wxf> <outDir>`, never a generated script. `Catch` + `Exit[1]` lives once in
> `driver.wls` (#561, #572: PSALTer throws uncaught and has a bare `Quit[]`), and the verdict is
> the artifact plus its sentinel line, never the exit status. It self-guards the lane with
> `pgrep` (the lane hook cannot see a Python-launched kernel) and tests the engine as a file
> (#565). A `wolframclient` session was measured and rejected: `DefField` never returns inside
> one (`interfaces_decision.md` §2.5, reproduced by the orchestrator 2026-09-15).

`tidalcosmo/derive/wolfram_driver.py`: engine-idle guard before every launch
(`pgrep -f -i 'wolframscript|WolframKernel|MathKernel'` — the single-license rule is
project-wide, not per-session), per-run log capture, wall-clock timing, configurable
timeout, strictly serial. Ported in shape from `tidal/cli/_derive.py::_run_wolframscript`
(line ~6438). No parallelism anywhere in this subsystem, ever.

> **⚠ Amendment (I-REM, 2026-09-09) — two measured requirements this driver inherits.**
> Both were established by I-526 (`stage1_measurements.md` §2.3, §2.5, §10.2) and routed to
> the orchestrator without landing anywhere. §1 of this document already points here
> ("durable coverage belongs in `tidalcosmo/derive/wolfram_driver.py` (§4.4), beside the
> engine-idle guard"); this is that coverage, stated.
>
> 1. **Never judge a `wolframscript` run by its exit status.** Wolfram 14.3 segfaults on
>    shutdown intermittently — *after* the work completes and the results are correct. The
>    CTEG Tier-1 gate run exited **143** having produced correct artifacts, and a trivial
>    one-line command produced a segmentation fault once in three runs, after printing the right answer.
>    **Decide success from what was produced**: the expected file exists, the expected
>    sentinel line is present. An exact capture of combined stdout+stderr is not reliable
>    either. (This is also what #541 fixed in `verify-wolfram-setup.sh`.)
>
> 2. **Set `QT_QPA_PLATFORM=offscreen` on every launch**, beside the engine-idle guard.
>    PSALTer exports a PDF through the Wolfram front end **unconditionally and without a time
>    limit**, and `$NoExport` does not guard it. When no Qt platform plugin can be
>    initialized, Qt aborts the front-end process and the export **blocks indefinitely rather
>    than failing** — measured at 20.001 s against a 20 s cap, twice; with the variable set,
>    3.8 s. **Exposure begins at `DefField`**, seconds into any run, not at
>    `ParticleSpectrum`. A hang is the dangerous failure here: it consumes the single-license
>    Wolfram lane indefinitely and is indistinguishable from "PSALTer is slow on this
>    theory" — which is the very measurement the Wave-2 cost run exists to make. The
>    variable belongs at this chokepoint deliberately, not container-wide.

## 5. The exporter

**`tidalcosmo/derive/wolfram/Stage1Export.wl`**, with a provenance header recording that
the structure-harvesting pattern is modeled on `SupplementalMaterials-2607`
`WolframLanguage/ParticleSpectroscopy/JuliaExport.m`, and that the coefficient-tensor and
explicit-label export exist nowhere upstream.

> **⚠ Amendment (I-526 live probe, 2026-09-07 — PSALTer `bb45adb0`; restated
> I-REM, 2026-09-09). The six-item list that stood here was WRONG**, omitting
> `$LocalWaveOperator` and `$LocalPropagator` on two implicit justifications that
> both fail on the live install:
>
> | test | result |
> |---|---|
> | `assoc[WaveOperator] === $LocalWaveOperator` | **False** — not redundant; omitting it discards information |
> | `$LocalPropagator` populated under `ShowPropagator -> False` | **True** (35 leaves) — not display-only, and not suppressed |
>
> Measured with the propagator display both on and off, identically. Evidence and
> transcript: `stage1_measurements.md` §7; tracked on #523/#527.
>
> **The fix is not a longer list.** A copied enumeration is what went wrong here, and
> a second copied enumeration would go wrong the same way at the next revision. The
> instruction below is therefore stated as a rule with a source of truth, and the
> count is recorded only as a dated measurement.

**Input surface** (§0.4): the theory association (`<Name>@WaveOperator`,
`@PseudoDeterminant`) plus **every** `xAct`PSALTer`Private`` global that
`ParticleSpectrum` populates en route.

**Enumerate them from the source, not from this document.** They are declared together at
`ParticleSpectrum.m:74-81` in the installed tree; read that at the revision you are pinned
to, and guard the resulting set with a fixture test that fails loudly if a name moves or one
is added — §0.4 already requires exactly that, because the exporter reads private symbols and
is therefore pinned to a PSALTer commit. **Measured 2026-09-09 on `bb45adb0`: eight**, which
§0.4 enumerates. Expect that number to change with the revision; do not hard-code it.

Read defensively — normalize keys, treat a bare unevaluated symbol as *absent*, and handle
both the populated case and the degenerate one (`UnitarityConditions =
Text["(Demonstrably impossible)"]` with an empty constraint matrix is a legal output).

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

   > **⚠ Amendment (I-REM, 2026-09-09 — PSALTer v2.0.2 `bb45adb0`, Wolfram 14.3).
   > `EXPIRES-WITH: #543` → EXPIRED 2026-09-11 (I-543).** Calibrating against the published
   > values and the committed upstream `.wxf` remains the right practice; the prohibition on
   > this install's output is lifted: with the two Function Repository dependencies registered
   > (`scripts/psalter/register_resources.wl`), the Vector theory's result association carries
   > both published blocks in both keys (`repro_543.wl B`: `½[Def²(Θ₂−Θ₁)−Θ₃]`, `½[−Def²Θ₁−Θ₃]`),
   > on 14.3.0 and on 14.2.1. What follows describes the pre-resolution state and is kept for
   > the record.
   >
   > ~~Calibrate against the published values and the committed
   > upstream `.wxf`; do NOT calibrate against what this install returns.~~
   >
   > On this install the Vector theory's **result association** — the object the reader
   > actually consumes — comes back as `{{-Def²Θ₁/2}}`: the `Θ₃` mass term and the `Θ₂`
   > contribution are both **absent**, while the rendered spectrograph showed the published
   > blocks. Measured on the published `VectorTheory.m`, identically with
   > `ShowPropagator` on and off; **#542**.
   >
   > So the gate as originally written would fail against a *correct* install, or — worse —
   > be quietly relaxed until it passed against degraded values, baking the defect into the
   > contract. The oracles to gate on are the author's own artifacts: the published
   > expressions, and `tests_cosmo/fixtures/psalter/ParticleSpectrographVectorTheory.wxf`,
   > which exists precisely so the reader can be tested without Wolfram installed
   > (`PROVENANCE.md`).
   >
   > **Open question this hands to #527, not for us to guess:** *which artifact carries the
   > published values* — the association keys, or what `ConstructSpectrograph` renders. They
   > are different objects and they disagree here. Settle it before specifying a gate
   > against either.
   >
   > **Do not generalize from "the render looked right".** That contrast holds for this
   > small probe only. The CTEG render on the same install is itself visibly broken —
   > unevaluated `Grid[...]` source, `$Failed`, hundreds of `Indeterminate`, ending in
   > "(Demonstrably impossible)"; artifact and counts at
   > `docs/cosmology/evidence/tier1-20260907/`. A plausible-looking PDF is not evidence of a
   > working install.
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

> **⚠ Amendment (I-526, 2026-09-07 — PSALTer `bb45adb0`).** The timestamped
> checkpoints described here **do not need to be built.** PSALTer already emits a
> stage trace: every function defined through its `StackSetDelayed` wrapper shells
> out an `echo` naming itself on entry whenever there is no notebook front end
> (`ConstructSpectrograph/CLICallStack.m:8-11`). A plain `wolframscript -file` run
> therefore prints one line per function entry, covering all seven `Construct*`
> stages. The Tier-1 gate timestamps that stream
> (`scripts/psalter/stamp_lines.py`) and extracts stage boundaries afterwards, so
> the run stays byte-for-byte unmodified while still producing a trace — 18,436
> stage events across 138 functions on CTEG. See `stage1_measurements.md` §4.2.
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
| oracle T1 | our `CTEG` run's association matches the committed `.mx` key-by-key. **PASSES since 2026-09-11** (`evidence/tier1-20260911-pass/`, `VERDICT: MATCH`; it failed until the two undeclared Function Repository resources were registered — #543), and a pass certifies only those two keys — see §3 |
| oracle T2/T3 | gauge-generator count exact; massive content with dictionary-exact `m²`; massless polarization count; `Reduce`-verified equivalence of unitarity conditions **including the massless one** |
| reader | `A23` blocks are dims (2,4,2) with `2·1+4·3+2·5 = 24`; `Vector` blocks equal the **published** expressions and the **committed upstream `.wxf`** exactly (symbolic difference simplifies to zero) — published values and the `.wxf` remain the *reference*; the 2026-09-09 prohibition on this install's own output is **lifted since certification** (2026-09-11, see the amendment below), and one open question survives: which artifact carries the published values — the association keys, or the `ConstructSpectrograph` render (§8, handed to I-S1A); placeholder/plural-key/degenerate cases unit-tested |
| reject rule | a bare-numeric term errors with a hint naming the term; the legacy `theory.toml` Lagrangian is rejected |
| exporter | in-kernel reconstruction `SameQ`; Python-side numeric agreement ≤ 1e-12 at 5 random rational coupling points; label calibration passes on ≥ 3 single-operator probes; `A[0]` vanishes **iff non-sampled constants are excluded from the coupling vector** — see the amendment below (coupling-linearity) |
| package + driver *(was "emitted script"; D-A, 2026-09-15)* | the manifest records signature and `ε` and every consumer refuses any other; the driver is fixed, never generated; no repo-absolute paths |
| cost | measured, not guessed, with checkpoints; one explicit go/no-go sentence; `Method` inertness confirmed live |

**Where tests live.** `scripts/run_wolfram_tests.sh` auto-discovers `tests/wolfram/test_*.wls`
by glob, so a new suite is picked up automatically — keep it **fast** (Vector-scale only,
seconds to minutes) or it will slow every future full-suite run. The oracle and PGT+EM runs
belong in `scripts/psalter/` wrappers, outside the glob. Python tests skip cleanly when
`wolframclient` is absent.

## 8. Risks

| risk | mitigation |
| --- | --- |
| ~~`UsingFrontEnd@Export` fails headless~~ **it HANGS headless, from `DefField` onward** | **Resolved 2026-09-07 (I-526):** set `QT_QPA_PLATFORM=offscreen`; the export then completes in ~4 s. Cause is that no Qt *platform plugin* could be initialized — `offscreen` was the only one with all its libraries present — **not** that `DISPLAY` is set, which was the first guess and is wrong. Exposure begins at `DefField`, not at the spectrograph. **Do not guard the export**; the guard would delete a working capability, and a hang (unlike a failure) consumes the single-license lane indefinitely. See the §0.6 amendments and `stage1_measurements.md` §2.3 |
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
   + driver**, text-emission only, with golden-text tests (§4) *(amended 2026-09-15, D-A: no
   generator and no golden-text tests — the package is unit-tested under `tests/wolfram/`)* → commit
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
   wrong in the other direction. Applied 2026-09-04, and extended by the coherence pass and
   the scientific review — **the live count is the ledger at the head of
   `spectrum_design.md`, not a number repeated here** (it said "six" and had already drifted).
   The notes cover
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
