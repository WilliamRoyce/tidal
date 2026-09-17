# I-S1B — Stage-1 Wolfram interface, exporter, cost run, and the `A[0]` answer

> **STATUS: HELD (outline) — Wave 2. Depends on D-C and I-S1A-core** (D-A, D-B and the
> conventions rule recorded 2026-09-15 and written in below). Written 2026-09-13 at Wave-1
> approval so nothing decided in Wave 1 is forgotten; the orchestrator updates the marked
> sections at each decision. **Not dispatched in Wave 1.** It becomes a
> full prompt, then READY, at Wave-2 planning.

| | |
|---|---|
| **Issue** | **#495** (WS6 umbrella) · #527 (I-S1A, the half that waited) · #522, #523 · #561 (the Wolfram-side `Exit[1]` deferred there) · #572 (PSALTer's bare `Quit[]`) · #566 (R-1) · #488 |
| **Milestone** | M-parallel |
| **Wave** | 2 |
| **Wolfram lane** | **Yes — the cost run and Tier 2.** One `wolframscript` at a time, machine-wide; `ensure_registered.sh` and `verify --require-psalter` exit 0 before and after. |
| **Depends on** | D-A, D-B ✅ (R-1, 2026-09-15); D-C (R-C); I-S1A-core merged (loader, store, launcher, dataclasses) |
| **Owned paths** | `tidalcosmo/derive/wolfram/` (the committed package, its exporter and the fixed `driver.wls`) · `tidalcosmo/cli/` — the `tidalcosmo derive` subcommand only · `tests/wolfram/` additions for it · `tests_cosmo/test_derive_wolfram*.py` · the cost-run record under `docs/cosmology/` |
| **NOT owned** | `tidalcosmo/{config,spectrum}/` beyond what the interface must call (I-S1A-core's; amend at the instruction site and report) · `tidalcosmo/{background,spectator,validity}/` · `pyproject.toml` |

### Decision dependencies (the orchestrator updates these before the Wave-2 rewrite)

| open decision | sections affected | what changes |
|---|---|---|
| **D-A** (R-1) — **recorded 2026-09-15: one committed package + a fixed driver** | the interface; the exporter; `Catch` + `Exit[1]`; the derivation's identity | The exporter is a **function in the committed package**, called by the **fixed `driver.wls`** (argv `<theory.wxf> <outDir> [--validate-only]`; sentinels `STAGE1 PASSED` / `STAGE1_ERROR=<tag>: <message>`; exit 0/1/2). Start from the R-1 prototype `scripts/research/interfaces/wolfram/{Stage1Proto.wl,driver.wls}` (it reproduced `vector_smoke.wls`'s wave operator `SameQ`, `interfaces_decision.md` §2.3). User text is parsed held and whitelisted before PSALTer runs (§2.4 — keep all eleven failure fixtures as `tests/wolfram/` cases); field and coupling symbols are created in `Global` (PSALTer derives contexts from `ToString`); the theory name is fingerprint-derived; `Catch` + `Exit[1]` lives once in the driver (#561, #572). The derivation's identity is the loader's **fingerprint** — there is no script to hash. |
| **D-B** (R-1) — **recorded: Option A′** | the input the package consumes; the CLI | The package receives the loader's canonical theory as WXF; nothing here reads a user file. **`tidalcosmo derive <run-or-theory>.yaml`** wires I-S1A-core's loader, launcher and store with `--test` (never starts a kernel; mirrors `cobaya-install --test`), `--force`, `--dry-run`, `--timeout` (default none) and `--validate-only`, as run in R-1's `derive_proto.py` (§3.3). **Schema gaps this prompt owns** (§3.8, each added under the schema-evolution rule): expansion geometry (which formulation supplies R and T from which fields), derived fields (`F = dA`), gauge-symmetry declarations (provisional), and widening the operator whitelist beyond the Vector theory's `CD`, fields and indices. |
| **D-C** (R-C) — **recorded 2026-09-17: O1′** | the shared toolbox; the post-Riemannian rewrite | The **rewrite is one committed function shared by both branches** — `ChangeCurvature[L, CDT, CD]` then the connection-difference tensor substituted by the contortion of the torsion, using the **corrected** identity (`−½(T^a{}_{bc} + T_b{}^a{}_c + T_c{}^a{}_b)`), **never legacy's text** (#582), carrying its three guards (antisymmetric part = torsion, metric compatibility, the Einstein–Cartan vector-torsion check) and gated against xMAG as an oracle. The toolbox stays PSALTer-shaped for this branch; the FRW branch adds xPand at M3. Curvature signs are xTensor's defaults, asserted per kernel (`conventions.md` §2.1). |
| **Conventions** (R-1) — **recorded: `conventions.md` canonical** | the exporter's `conventions` block | The user writes in PSALTer's convention, so **nothing is mapped**: the spectrum records `{"signature": [1,-1,-1,-1], "epsilon0123": 1}` (`PSALTer.m:101`, `DefGeometry.m:69-72`), the manifest repeats it and every reader refuses others. The spectrum-side check is a reference theory with a published healthy/ghost verdict run through the whole pipeline (`conventions.md` §5). |
| **I-S1A-core** | everything that consumes the dataclasses, the driver and the validator | Import paths from I-S1A-core's report. *(To be filled.)* |

## Outline — what this prompt will ask for

1. **The committed Stage-1 package and the fixed `driver.wls`** (D-A), launched by I-S1A-core's
   `tidalcosmo/derive/launcher.py`; the artifact verdict, the lane self-guard and the file-based
   engine test are already in the launcher — do not duplicate them. Plus the **`tidalcosmo derive`
   CLI** (D-B).
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
6. **#561's Wolfram-side half**: `Catch` at the top level of `driver.wls` with `Exit[1]` on an
   uncaught `Throw`; PSALTer's bare `Quit[]` (#572) is caught by the launcher's artifact check; a
   test that watches each fail first.
7. **The schema gaps named in the D-B row**, each added as an optional key under the
   schema-evolution rule (register row *Theory loader is permanent*), with a torsion theory run
   through the package end to end.

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
