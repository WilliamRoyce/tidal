# `docs/cosmology/handoffs/` — the prompt index

Every session that did work for the cosmology programme (#488) was dispatched with a prompt
file kept here. **The prompts are kept, not deleted**, because a prompt records *what was
asked*, which is what lets a later reader judge whether the deliverable answered it.

**None of these is an open assignment**, whatever their imperative mood — H2 in particular
opens "if it does not exist yet, ask before proceeding", which was live instruction text for a
task finished five days earlier. Read the status column here first.

**This index is not the state of the programme.** The **wave board** in
[`../../COSMOLOGY_PROGRAM.md`](../../COSMOLOGY_PROGRAM.md) is the single source of truth for
what is merged, dispatched or planned; this file only maps prompts to their outputs.

## Two series

- **H-series — research.** Executed during planning, before implementation began. All
  complete. Their prompts carry a free-text `**Program:** / **Mode:** / **Artifact:**` header
  and a `> **STATUS: COMPLETE**` blockquote, with no wave or lane field (waves did not exist
  yet).
- **I-series — implementation.** Dispatched per issue under the delegation protocol. Each
  carries a uniform header table: issue, wave, **Wolfram lane**, dependencies, owned paths.
- **R-series — wave research** (from Wave 1). Research memos dispatched *inside* a wave under
  the same protocol and header table as the I-series, because the builds behind them are held
  on their decisions. `R-1` decides D-A/D-B; `R-C` decides D-C.

**Statuses, from Wave 1 on** (dispatch-cadence rule, 2026-09-13): every prompt of a wave is
written at approval. `> **STATUS: HELD — depends on …**` means written, waiting on a decision
named in its dependency table — **do not dispatch**. `> **STATUS: READY**` means every
dependency is recorded and the text updated — the user dispatches. Wave-2 outlines are held
prompts with named gaps, rewritten in full at Wave-2 planning.

**H7 has no prompt file.** It was executed inside the orchestrator's own planning session
rather than dispatched, so only its output was committed. That is why this directory holds
seven `H*.md` files while the programme document speaks of eight handoffs.

## H-series — research (all complete)

| Prompt | Task | Output |
| --- | --- | --- |
| `H1.md` | TorC pipeline audit — settled O1's scope and how the CAMB patch is made (#498) | `../torc_pipeline_audit.md` |
| `H2.md` | Observable-ladder feasibility — established that O2 and O3 are *different numerical problems* (#500–#510) | `../observable_ladder.md`, `../magnetic_field_background.md` |
| `H3.md` | Solver design study — two engines over one shared core; matrix-WKB designed and prototyped (#517–#520) | `../solver_design.md` |
| `H4.md` | New-package design — the strangler-fig migration and the `tidalcosmo/` scaffold (#513–#516) | `../repo_reshape.md` |
| `H5.md` | Literature acquisition — 20/20 fetched and title-verified (#497) | `literature/`, `../../references.md` |
| `H6.md` | Numerical polology design — the two-stage spectrum architecture and its primary algorithm | `../spectrum_design.md` |
| *(H7)* | **No prompt file.** Spectator-route scope, executed during planning | `../spectator_route.md` |
| `H8.md` | Stage-1 engineering **study** — six live-source findings correcting H6 (#521–#523) | `../stage1_engineering_plan.md`, `scripts/research/psalter_stage1/` |

## I-series — implementation

| Prompt | Issue(s) | Wave | Wolfram lane | Status | Output |
| --- | --- | --- | --- | --- | --- |
| `I-524.md` | #524 (M0) | 0 | no | **merged** | `tidalcosmo/` installable, extras, console script, CI lane |
| `I-525.md` | #525 (M0.5) | 0 | no | **merged** | `scripts/oracles/`, 185 frozen fixtures under `tests_cosmo/data/oracles/`; filed #535–#538 |
| `I-526.md` | #526 | 0 | **yes** | **merged**; gate resolved by I-543 (2026-09-11) | PSALTer v2.0.2 installed; Tier-1 MISMATCH on merge → #543. `../stage1_measurements.md`, `../evidence/tier1-20260907/` |
| `I-REM.md` | #545, #546, #540 | 0-completion | no | **merged** 2026-09-09 (#550) | instruction-site amendments, docs index, `oracle.yml`, tooling |
| `I-533.md` | #533 | 0-completion | no | **merged** 2026-09-11 (#552) | retire the M0 drop rows (`sweep`, `sample`, `analyze`, `plot`) |
| `I-543.md` | #543 (+#542, #549) | 0-completion | **yes** | **merged** 2026-09-11 (#557) — **gate passes, install certified** | resolve the Tier-1 gate — pass it, or locate the mechanism |
| `I-ONB.md` | #559 (#488, #543 hardening) | 0-completion | no | **merged** 2026-09-12 (#564) | one onboarding path: a fresh user with their own Wolfram ID lands in the certified configuration; fixes the fresh-host container failure |
| `R-1.md` | #566 (M1a; decides for #532, #527) | 1 | **yes** (≤ 30 min) | **merged** 2026-09-15 (#568) — D-A, D-B, conventions and a permanent loader decided with the user | `../interfaces_decision.md` (§3.8 schema stability), `../conventions.md` (canonical), `../r1_planning_record.md`; prototypes under `scripts/research/interfaces/`, lane evidence re-run by the orchestrator; issues #569–#573 |
| `R-C.md` | #567 (M3; informs #500, #501, #504) | 1 | **yes** (≤ 2 days + follow-up) | **reported** 2026-09-16 (PR #574, draft) — recommends O1′; follow-up lane on xMAG/xPand in progress | `../perturbation_tooling.md`, `../rc_planning_record.md`; probes under `scripts/research/perturbations/`; issues #575–#586 |
| `I-532.md` | #532 (M1a) | 1 | no | **HELD** — waits for R-C (its R-1 inputs recorded 2026-09-15) | `tidalcosmo/{background,spectator,validity}/`: protocol, CAMB seam, pass-through, flag schema, `../veto_placement.md` (group practice first), reference oracles + convention record, `camb` core |
| `I-S1A.md` | #527 (M-parallel; #522, #523) | 1 | no | **HELD** — I-532's flag type (rewritten 2026-09-15 for D-A, D-B and a permanent loader) | `tidalcosmo/{config,derive,spectrum}/` Python only: the permanent Option A′ loader, the spectra store, the driver launcher, WXF reader + Stage-2 dataclasses; `A[0]` as a named unknown |
| `I-S1B.md` | #495 (M-parallel) | 2 | **yes** | **HELD (outline)** — D-C, I-S1A-core | the committed Stage-1 package + fixed `driver.wls`, exporter, `tidalcosmo derive` CLI, cost run, the `A[0]` answer, Tier 2, #561/#572, the schema gaps of R-1 §3.8 |
| `M1b.md` | #494, #515 (M1b) | 2 | no | **HELD (outline)** — I-532, I-S1A-core | `SpectatorTheory(Theory)` on the permanent loader and store (`get_modified_defaults`, `ComponentNotInstalledError`, fingerprint as version), the `__init__` re-export, veto wiring per I-532's memo, ΛCDM identity gate, #515 |

## The Wolfram lane

**One `wolframscript` machine-wide**, so at most one lane-flagged session may run at a time —
the orchestrator included. Only a prompt whose header says **Wolfram lane: yes** may start a
kernel; the rest must not, and from Wave 1 each no-lane prompt carries the *complete list* of
commands that would start one (a fence that names a tool cannot see a wrapper script that runs
it, #555). Lane occupancy in Wave 1 is sequential by construction: R-1 briefly, then R-C.

## Conventions worth knowing before reading one

- A delegate works in its **own git worktree** off `feat/cosmology-program`, branches
  `cosmo/i<issue>-<slug>`, **never merges**, and opens a **draft PR at its first commit** —
  `test.yml` fires on `pull_request`, so a branch with no PR is never seen by CI.
- **Delegates never version-bump, tag, or edit the changelog.** The orchestrator bumps once
  per wave; concurrent bumps guarantee a `pyproject.toml` conflict and a tag collision.
- Each prompt carries a **scope fence** naming the adjacent temptations, and a **flaw
  protocol**: amend a wrong design *at the instruction site* and report it; **stop and report**
  if the design cannot be built as specified.

Full protocol, including the merge checklist and the wave-boundary checklist:
[`../../COSMOLOGY_PROGRAM.md`](../../COSMOLOGY_PROGRAM.md) §"Implementation delegation
protocol".
