# I-S1B — Stage-1 Wolfram interface, exporter, cost run, and the `A[0]` answer

> **STATUS: HELD (outline) — Wave 2. Depends on D-A, D-B, D-C and I-S1A-core.** Written
> 2026-09-13 at Wave-1 approval so nothing decided in Wave 1 is forgotten; the orchestrator
> updates the marked sections at each decision. **Not dispatched in Wave 1.** It becomes a
> full prompt, then READY, at Wave-2 planning.

| | |
|---|---|
| **Issue** | **#495** (WS6 umbrella) · #527 (I-S1A, the half that waited) · #522, #523 · #561 (the Wolfram-side `Exit[1]` deferred there) · #488 |
| **Milestone** | M-parallel |
| **Wave** | 2 |
| **Wolfram lane** | **Yes — the cost run and Tier 2.** One `wolframscript` at a time, machine-wide; `ensure_registered.sh` and `verify --require-psalter` exit 0 before and after. |
| **Depends on** | D-A, D-B (R-1); D-C (R-C); I-S1A-core merged (dataclasses, driver, validator) |
| **Owned paths** | `tidalcosmo/derive/wolfram/` (the interface and the exporter, in whatever shape D-A chose) · `tests/wolfram/` additions for it · `tests_cosmo/test_derive_wolfram*.py` · the cost-run record under `docs/cosmology/` |
| **NOT owned** | `tidalcosmo/{config,spectrum}/` beyond what the interface must call (I-S1A-core's; amend at the instruction site and report) · `tidalcosmo/{background,spectator,validity}/` · `pyproject.toml` |

### Decision dependencies (the orchestrator updates these before the Wave-2 rewrite)

| open decision | sections affected | what changes |
|---|---|---|
| **D-A** (R-1) | the interface's shape; the exporter's packaging; where #561's Wolfram-side `Catch` + `Exit[1]` lives; how `derivation_hash` is computed | Generated `.wls` → a templater in Python with golden-text tests and the exporter as a committed `.wl` the script `Get`s; a data-taking package → the exporter *is* the package entry, `tests/wolfram/` unit tests, a three-line data-loading stub; a `wolframclient` session → the package plus session-lifecycle code in the driver, and error signalling through the session rather than exit status. *(To be filled: which, and R-1's prototype tables.)* |
| **D-B** (R-1) | the input model the interface consumes | The interface reads I-S1A-core's dataclasses, which load D-B's format; nothing here reads a user file directly. *(To be filled: any derivation-only fields.)* |
| **D-C** (R-C) | what the shared convention-free toolbox must carry for the FRW branch | If R-C recommends building on an xAct-family package, the toolbox's field declarations and term structure must be emittable in that package's form as well as PSALTer's; if build-own, the toolbox stays PSALTer-shaped and the FRW branch adds its own emitter at M3. *(To be filled.)* |
| **Carriage rule** (R-1) | the exporter's `conventions` block; the mapping applied before PSALTer sees the Lagrangian | *(To be filled: the user-facing convention, PSALTer's signature and `ε₀₁₂₃` from source, the mapping's test.)* |
| **I-S1A-core** | everything that consumes the dataclasses, the driver and the validator | Import paths from I-S1A-core's report. *(To be filled.)* |

## Outline — what this prompt will ask for

1. **The Stage-1 Wolfram interface in D-A's shape**, driven by I-S1A-core's
   `wolfram_driver.py`; the never-judge-by-exit-status rule and the file-based engine test
   are already in the driver — do not duplicate them.
2. **The exporter** (`stage1_engineering_plan.md:565-570`, its packaging per D-A): emits the
   Stage-2 contract the I-S1A-core dataclasses read — schema version, theory, PSALTer commit,
   `conventions`, ordered couplings, `J^P` labels, `A[(n+1) × 3 × dim × dim]`, gauge rank,
   partition, source constraints, gauge generators, polarizations, massive spectrum, the
   four-outcome unitarity status, `z_degree`, `k4_guard`. The `conventions` block records
   PSALTer's signature and `ε₀₁₂₃` as R-1 found them.
3. **The `A[0]` answer.** Whether non-sampled constants enter the coupling vector — decided
   from a real export on the Vector theory and one theory with a fixed constant; the
   dataclasses' named unknown is resolved and its docstring amended at the instruction site.
4. **The cost run** (`stage1_engineering_plan.md` §0.5): wall time and peak memory for the
   Vector theory and for `A23` through the interface, recorded with the certified
   configuration's fingerprint; the number the sampling-time budget stands on.
5. **Tier 2** of the PSALTer gate on the interface's own output — the same fixtures I-S1A-core
   asserts against, produced by this path, byte-compared.
6. **#561's Wolfram-side half**: `Catch` at the top level with `Exit[1]` on an uncaught
   `Throw`, placed where D-A's shape puts it; a test that watches it fail first.

## Success criteria — to be finalized at the Wave-2 rewrite

- The interface produces, for the Vector theory, a contract whose `WaveOperator`-derived
  blocks equal I-S1A-core's committed fixture to the digest.
- `A23` through this path yields `(2,4,2)` and `24`; digest asserted.
- `A[0]` is answered with the export that decided it committed as a fixture.
- Cost-run numbers recorded with the fingerprint; lane time within the recorded budget.
- The `Exit[1]` test fails before the fix and passes after; `tidal derive`'s artifact check
  (#561, Python side) still covers the exit-0-no-output case.
- `verify --require-psalter` exit 0 after; `CI <run-id>: success`; pyright; boundary test.

## Scope fence — provisional

No Stage-2 evaluator. No FRW branch (M3). No change to the PSALTer pin. Nothing under
`tidalcosmo/{background,spectator,validity}/`. No `pyproject.toml` edits.

## Working rules

The standard set from `docs/COSMOLOGY_PROGRAM.md` §"Prompt template" — worktree off
`feat/cosmology-program` (`git worktree add /tmp/tidal-is1b -b cosmo/is1b-stage1-wolfram
feat/cosmology-program`); draft PR at the first commit; never merge, bump, tag or edit the
changelog; assertions verified before they land; a verification never writes into what it
verifies; guards carry `EXPIRES-WITH:`; one kernel machine-wide; no environment-specific
absolute paths; conventional commits; American English; no attribution trailers.

## Report back

Branch · PR · `CI <run-id>: <conclusion>` · each criterion with its output · the `A[0]`
answer and the fixture that decided it · the cost-run table · both verify outputs ·
amendments made · discoveries to route.
