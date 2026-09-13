# I-S1A-core — Stage-1 Python side, minus the Wolfram interface

> **STATUS: HELD — depends on D-A, D-B, the convention-carriage rule (all from R-1) and on
> I-532 having merged (its flag type).** Written 2026-09-13 at Wave-1 approval; the
> orchestrator updates the marked sections at each decision and marks this READY after I-532
> merges. **Do not dispatch while HELD.**

| | |
|---|---|
| **Issue** | **#527** (I-S1A) · #495 (WS6 umbrella) · #522 (the validator's origin), #523 (the harvest surface) · #488 |
| **Milestone** | M-parallel |
| **Wave** | 1 — build, second; **after I-532 has merged** |
| **Wolfram lane** | **NO.** Never start a kernel — every test runs against the committed `.wxf` fixtures. The complete list of things that would start one: `wolframscript`, `tidal derive`, every `examples/*/run.sh`, `scripts/psalter/{run_tier1_gate.sh,ensure_registered.sh,register_resources.wl,repro_543.wl,vector_smoke.wls,probe_*.wls}`, `scripts/verify-wolfram-setup.sh`, `scripts/run_wolfram_tests.sh`, `scripts/full_test.sh`, the `/derive` skill (#555). |
| **Depends on** | D-A, D-B, the carriage rule recorded; I-532 merged and its flag type's import path published |
| **Owned paths** | `tidalcosmo/{config,derive,spectrum}/` — **Python only, never `tidalcosmo/derive/wolfram/`** · `tests_cosmo/test_{config,derive,spectrum}*.py` · `tests_cosmo/fixtures/psalter/` (read-only) |
| **NOT owned** | `pyproject.toml` (I-532 owns it this wave; if D-A chose `wolframclient` and it is absent, the orchestrator adds it before you start) · `tidalcosmo/{background,spectator,validity}/` · any `.wl` or `.wls` · the design docs |

### Decision dependencies (the orchestrator updates these before marking READY)

| open decision | sections affected | what changes |
|---|---|---|
| **D-A** (R-1) | Deliverable 2, the driver | Generated `.wls` → the driver runs a script file; a data-taking package → the driver invokes a committed entry with a data file; a `wolframclient` session → the driver owns the session lifecycle. Guard, offscreen, artifact-verdict and refusal are the same in all three. *(To be filled: which.)* |
| **D-B** (R-1) | Deliverable 1, the loader | Loads D-B's format directly. *(To be filled: the block's fields, and whether the derivation reads a sub-block of the Cobaya YAML.)* |
| **Carriage rule** (R-1) | Deliverables 1 and 3 | The convention the user writes in; the mapping into each backend; what `conventions.signature` the fixtures carry. *(To be filled.)* |
| **I-532's flag type** | Deliverable 3 | The verdict uses it; import path from I-532's report. *(To be filled.)* |

## Why this exists

The two-stage spectrum architecture (`spectrum_design.md` §1; `stage1_engineering_plan.md`)
runs Wolfram/PSALTer once at derivation time and evaluates the exported contract numerically
per sample, with **no Wolfram at sampling time**. Its Python side has four pieces that do not
depend on how the Wolfram interface is shaped, and they are yours:

1. **The input model and its coupling-linearity validator.** PSALTer does **not** enforce
   coupling-linearity — `ParticleSpectrum::NonLinearCouplings` is defined but never thrown, and
   a bare numeric coefficient passes silently (#522). `spectrum_design.md` §4.1's "reject, never
   auto-assign" rule is therefore entirely ours and **load-bearing for correctness**.
2. **The `wolframscript` driver**, which carries three measured requirements (#543,
   `stage1_engineering_plan.md:534-563`): one kernel machine-wide; `QT_QPA_PLATFORM=offscreen`
   on every launch (without a loadable Qt platform plugin `UsingFrontEnd@Export` hangs forever,
   from `DefField` on); and **never judge a run by its exit status** — 14.3 segfaults on
   shutdown after producing correct results (the certifying gate run exited 143), and an
   uncaught `Throw` exits 0 (#561). Decide from the artifact.
3. **The WXF reader and the Stage-2 contract dataclasses**, tested against the committed
   `.wxf` answer keys (`tests_cosmo/fixtures/psalter/`, `PROVENANCE.md`).
4. **`A[0]` — represented, not answered.** Whether non-sampled constants enter the coupling
   vector needs a real export (I-S1B). The dataclasses represent both forms and carry it as a
   named unknown.

What is **not** yours: the Wolfram interface itself (D-A's shape) and the exporter `.wl` —
those are I-S1B's, built together so the session that writes them also runs them on the lane.

## Ordered reading list

1. `docs/COSMOLOGY_PROGRAM.md` — register rows *I-S1A split*, *Convention carriage*, *Config
   layer* (provisional, disposable), *`A[0]` owner*, *Dispatch cadence*.
2. `docs/cosmology/interfaces_decision.md` (R-1) — D-A, D-B, the carriage rule.
3. `docs/cosmology/stage1_engineering_plan.md` §4.2 (`:481-500`), §4.4 (`:532-563`), §5
   (`:565-640`), §0.5 (`:177-216`), the gates at `:753-764`.
4. `docs/cosmology/spectrum_design.md` §4.1 (`:238-261`), §6.1, `:346-366`, `:610-613`.
5. `tidalcosmo/{config,derive,spectrum}/README.md` (sketches).
6. `tests_cosmo/fixtures/psalter/PROVENANCE.md` (which artifact carries the published values).
7. I-532's report (the flag type's import path).
8. The pinned PSALTer sources: `Sources/DefField.m:11-24` (symmetry classes),
   `ParticleSpectrum.m:74-81` (the private globals, **eight** on `bb45adb0` — enumerate,
   never hard-code).

## Deliverables

1. **`tidalcosmo/config/`** — frozen dataclasses (`Theory`, `Sector`, coupling vocabulary per
   `repo_reshape.md:494-497`) loading **D-B's format directly** *(to be filled)*. The
   Lagrangian is read in the **one documented convention** R-1 chose; the mapping into each
   backend's fixed convention is a pipeline rule with a test; every emitted artifact records the
   convention it is in. Every docstring: *provisional / disposable*, `EXPIRES-WITH: #490` —
   M1b **replaces** this layer, it does not extend it. Validator: field symmetry classes exactly
   PSALTer's registered list (`DefField.m:11-24`; torsion is
   `RegisterFieldRank3Antisymmetric23`), naming the legal values; coupling-per-operator reject
   rule with `error_with_hint`-style hints (the value belongs in the coupling, not the term);
   echo PSALTer's own `UnknownCoupling` text; the legacy
   `examples/torsion_gertsenshtein/theory.toml` Lagrangian *string* as the committed negative
   fixture — embed the string, not the path; assert no background quantity leaks into the
   vacuum block.
2. **`tidalcosmo/derive/wolfram_driver.py`** — the engine test is a **file**
   (`$HOME/.local/wolfram/engine/<series>/Executables/WolframKernel`; `command -v
   wolframscript` is satisfied by the image's cloud-only client, #565); `QT_QPA_PLATFORM=
   offscreen` at this chokepoint; success decided from artifacts, never exit status; refusal
   unless `verify-wolfram-setup.sh --require-psalter` exited 0 (cached per process); strictly
   serial; timeout; log capture. Shaped to D-A *(to be filled)*. Unit-tested with subprocess
   mocks — **no kernel is ever started by this prompt**.
3. **`tidalcosmo/spectrum/`** — the WXF reader (defensive key normalization; bare symbol =
   absent, `PROVENANCE.md:19-23`) and the Stage-2 contract dataclasses per
   `spectrum_design.md` §6.1 / `stage1_engineering_plan.md:605-620`: schema version, theory,
   PSALTer commit, `conventions` (asserted against the expected backend convention — a
   mismatch fails), ordered couplings, per-sector explicit `J^P` labels, the coupling-linear
   coefficient tensor `A[(n+1) × 3 × dim × dim]`, gauge rank, the massive/massless partition,
   source constraints, gauge generators, massless polarizations, massive spectrum, unitarity
   status as the **four-outcome** enum (closed form / provably impossible / absent / PSALTer
   halted), `z_degree`, `k4_guard`. `A[0]` representable in both forms with a named unknown and
   a test that exercises both. A test that reads `ParticleSpectrum.m:74-81` from the installed
   tree **when present** and compares the private-global set to the reader's — skipping cleanly
   otherwise. The verdict uses I-532's flag type *(import path to be filled)*.

## Success criteria — verified from artifacts

1. `A23` block dims `(2,4,2)` with `2·1 + 4·3 + 2·5 = 24`; `Vector` blocks equal the published
   expressions from the committed fixture (state **which artifact** carried them,
   `PROVENANCE.md:87-93`); fixture digests asserted.
2. A convention mismatch between a fixture and the expected backend is a failing test.
3. ≥ 1 test per validator message; the legacy Lagrangian string is rejected.
4. Driver tests cover: the file-based engine guard, offscreen in the environment, the
   artifact-based verdict (exit 0 with no artifact → failure; exit 143 with a good artifact →
   success), the verify-refusal.
5. `import tidalcosmo` ≤ 600 ms; pyright clean; `tests_cosmo/test_package_boundary.py` green.
6. Draft PR from the first commit; `CI <run-id>: success`; `ruff`, `ruff format --check`,
   `cspell`.

## Scope fence

No Wolfram interface of any shape — no generator, no golden-text tests, no `.wl`, no
`.wls`. No Stage-2 evaluator. No kernel. No new user-facing file format beyond D-B. No
`pyproject.toml` edits. The config layer is **disposable** — do not let it accrete a schema
that would constrain M1b's real one.

## Working rules

- Worktree off `feat/cosmology-program`: `git worktree add /tmp/tidal-is1a -b cosmo/is1a-stage1-python feat/cosmology-program`. **Never merge.**
- **Never version-bump, tag, or edit the changelog.**
- **Draft PR into `feat/cosmology-program` at your first commit**; CI is the gate; `CI <run-id>: <conclusion>`.
- Stay inside your owned paths.
- **Assertions get verified before they land; hypotheses do not have to be.**
- **A verification command must not write into what it verifies.**
- **Any guard you add carries its expiry** (`EXPIRES-WITH:`) — the private-globals guard expires with the PSALTer pin; the config layer with #490.
- **Read the tool's own known-issues before hypothesizing** (PSALTer's README known bug #1 is the one Wave 0 spent a week attributing to the engine).
- No kernel, ever (the list in the header).
- No environment-specific absolute paths in anything committed.
- Conventional commits; American English; no attribution trailers.

## If you find the design wrong

Amend a design-document error at the instruction site and report it. An architectural
contradiction: stop and report.

## Report back

Branch · PR · `CI <run-id>: <conclusion>` · each criterion with its output · where the `A[0]`
unknown is named (a docstring path, not only this report) · amendments made · discoveries to
route · suggested next step.
