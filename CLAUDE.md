# TIDAL: Tensor Integration and Derivation for Any Lagrangian

Symbolic physics pipeline: Lagrangian (xAct/Mathematica) -> JSON -> native PDE solver (SUNDIALS IDA/CVODE, leapfrog, scipy) with numpy spatial operators. All PDEs derive from the Lagrangian via symbolic computation -- never manually hardcode equations. Operates exclusively in the linearized regime (quadratic Lagrangians, linear PDEs).

## Project Structure

- `tidal/wolfram/` -- Wolfram pipeline modules (EulerLagrange.wl, ComponentDecompose.wl, ExportJSON.wl, CommonUtilities.wl, GaugeFix.wl)
- `tidal/solver/` -- PDE solver backends (ida.py, cvode.py, leapfrog.py, fields.py, operators.py, grid.py, coefficients.py, rhs.py, state.py, constraint_solve.py)
- `tidal/symbolic/` -- Python symbolic pipeline (_derive.py, json_loader.py)
- `tidal/cli/` -- CLI entry points (11 subcommands: derive, simulate, measure, inspect, list, validate, plot, sweep, analyze, sample, doctor)
- `tidal/inference/` -- Bayesian inference (priors, likelihood, constraints, MC, nested sampling via dynesty/PolyChord)
- `tidal/measurement/` -- Physics measurements (energy, conversion, mixing, spectra)
- `examples/` -- 19 physics examples (1+1D through 3+1D), each with theory.toml + .wls + data/*.json
- `research/` -- General quadratic PGT+EM Lagrangian enumeration (xAct/xTras scripts, TeX document, classification JSONs)
- `tests/` -- ~2,400 Python tests (76 files) + ~260 Wolfram test cases (tests/wolfram/*.wls)
- `docs/` -- Architecture and program docs (`docs/README.md` is the index; `docs/tex/` holds the technical reference)

## Key Commands

- `uv run pytest tests/ -x -q` -- Run Python tests
- `./scripts/full_test.sh` -- Full test suite (Python + Wolfram)
- `uv run tidal derive examples/<name>/theory.toml` -- Derive PDEs from Lagrangian
- `uv run tidal simulate examples/data/<name>.json` -- Run simulation
- `uv run tidal sweep examples/data/<name>.json --sweep "param=start:stop:N" --measure conversion --output sweep_out` -- Run parameter sweep
- `uv run tidal sample examples/data/<name>.json --prior "param=uniform:lo:hi" --likelihood "P_max:maximize" --method mc --n-samples 100 --output sample_out` -- Bayesian inference (MC or nested sampling)
- `uv run ruff check` / `uv run ruff format` -- Lint / format
- `uv run pyright` -- Type checking. **Requires `uv sync --all-extras` first.** Without the
  optional extras installed, pyright reports phantom unresolved-import errors for `jax` in
  `tidal/solver/modal_jax.py` -- CI syncs all extras and reports 0 errors, so a local-only
  failure there means the venv has drifted, not that the code is broken.

## Critical Conventions

- Pipeline: TOML config -> .wls generation -> wolframscript -> JSON spec -> Python solver
- Solver selection: IDA (DAE/constraints), CVODE (adaptive ODE), leapfrog (symplectic), scipy (general)
- Always check xAct symbol existence before defining: `If[!xTensorQ[M2], DefManifold[...]]`
- Parenthesize multiline Lagrangians in .wls files
- Use `DefConstantSymbol` for mass/coupling constants (not bare Symbol)
- **Constant names must not contain underscores** — Mathematica parses `X_Y` as `Pattern[X, Blank[Y]]`, corrupting symbolic computation. Use `mPhi2` not `m_phi_2`, `Bpeak` not `B_peak`.
- Cross-field decomposition requires passing `additionalFields` to `DecomposeToComponents`
- Background fields declared via `[[background_fields]]` TOML section
- Gauge fixing via `[[gauge]]` TOML section (presets: Lorenz, de Donder, Coulomb, temporal, axial)
- Velocity naming: v_{field_name} (e.g., v_phi_0, v_A_1) — E-L velocity form, not canonical momenta
- **User-facing errors must include hints**: Use `error_with_hint(msg, hints)` from `tidal.cli._console` instead of bare `error()` for all CLI error messages. Each hint should be an actionable suggestion (example syntax, available options, related commands, troubleshooting steps). See existing ~60 error sites across CLI modules for the pattern.
- **No environment-specific absolute paths in committed files**: TIDAL is meant to be distributed, so anything hardcoding an absolute container path under `/workspaces`, a user home under `/home`, a machine-specific venv, or a Claude project slug works only for its author. Derive the root instead — shell scripts use `REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"` (see `scripts/hpc_shuttle.sh:13`); devcontainer lifecycle commands run *from* the workspace folder, so plain relative paths are correct; Claude Code hooks use `$CLAUDE_PROJECT_DIR`; the Claude project-dir slug is the absolute path with every non-alphanumeric character replaced by `-`, so derive it with `sed 's/[^a-zA-Z0-9]/-/g'` rather than pasting it. Machine-local files (`.claude/settings.local.json`, gitignored) are the one exception.

## Workflow Rules

- **After completing any code change**, run relevant tests before moving on. Source→test mapping: `tidal/solver/X.py` → `tests/test_solver_X.py`, `tidal/cli/_X.py` → `tests/test_cli.py`, `tidal/measurement/` → `tests/test_measurement.py`. Unsure → full suite: `uv run pytest tests/ -x -q`
- **After completing a feature/fix**, commit promptly with conventional format (feat:/fix:/refactor:/test:/docs:). No Co-Authored-By trailer. Separate unrelated changes into distinct commits.
- **Fix lint/type/spell errors immediately** — `uv run ruff check --fix && uv run ruff format` after code changes. Fix pyright errors. Add domain terms to `cspell.json`, fix genuine typos.
- **Wolfram pipeline integrity**: ALL symbolic processing stays in Wolfram — never post-process equations in Python. Never skip/bypass the canonical pipeline; fix root causes.
- **Run long commands in background**: Use `run_in_background: true` on the Bash tool for any command that takes more than a few seconds — derivations (`tidal derive`), simulations (`tidal simulate`), sweeps (`tidal sweep`), and full test suites (`pytest tests/`). Continue other work while waiting; you'll be notified on completion. Do NOT poll or sleep.
- **Only ONE wolframscript at a time** — single engine license. NEVER run `tidal derive` in parallel.
- **Use minimal test theories** (scalar_field, coupled_scalars) before expensive derivations.
- **Negative energies** may be physical with (-,+,+,+) metric convention — don't "fix" without understanding the physics.
- **Before context compaction**, update all relevant docs and memory files.
- **Version bump after completing work**: After committing a completed feature/fix (all tests pass, no remaining tasks), bump the version: `--patch` for fixes/small changes, `--minor` for new features. Skip if mid-feature or WIP. NEVER bump the major version automatically. Use `python scripts/bump_version.py --{level} --commit --allow-dirty`. Default to bumping — only skip if you are about to make another commit immediately in the same sitting.
- **Update documentation after completing work**: After committing a feature/fix, identify and update relevant docs. Technical docs live in `docs/tex/*.tex` (the primary documentation). Project management (roadmaps, checklists) lives in `docs/*.md`. To find affected docs: `grep -rl "keyword" docs/tex/ docs/*.md`. Common update patterns: phase/issue status changes → `docs/ROADMAP.md`, `docs/NEXT_PHASES.md`; implementation substep done → active checklist in `docs/`; performance changed → whichever `docs/tex/*.tex` has benchmark tables for that subsystem; new error pattern → `docs/tex/troubleshooting.tex`; algorithm/architecture changed → whichever `docs/tex/*.tex` describes that component. See `docs/README.md` for the full documentation index. The `/sync-docs` skill checks for drift — run it after major features. Commit doc updates separately: `docs: update {topic} documentation`.
- **Only commit YOUR changes**: Before staging files, verify each file's diff contains changes YOU made in THIS session. Never stage files modified by parallel agent sessions or other worktrees. If `git status` shows unexpected modified files, check `git diff <file>` before staging.
- **Create GitHub issues proactively**: When you encounter bugs, improvement opportunities, technical debt, or notable discoveries during work, create a GitHub issue via `gh issue create` to build a searchable trail. This applies even for things you fix immediately — create the issue, then close it with the fix commit (`gh issue close N -c "Fixed in <commit>"`) so there's a record of what was found and how it was resolved. Always check for duplicates first: `gh issue list -S "keyword"`. Use appropriate labels from the existing set (bug, enhancement, documentation, validation, etc.). Include: clear title, context, relevant file paths, and why it matters. Skip only if: truly trivial (typo, formatting) or a duplicate already exists. **NEVER include any "Generated with Claude Code" footer or attribution in issue bodies.**

## Physics Coding Patterns

- **Specify success criteria before coding**: "Modal solver must agree with CVODE to RMS < 1%" — not just "implement modal solver". Include quantitative thresholds.
- **Wolfram derivations**: Read an existing .wls template first, generate new by modifying template, review diff against template before running wolframscript
- **After derivation**: Verify JSON has `canonical.hamiltonian_terms` — without this, all energy measurements fail silently. Run `tidal validate <json> --stability`.
- **Record derivation timing**: After a successful `tidal derive`, update the theory TOML header comment with `# Derivation timing: ~Xm wall (last verified: YYYY-MM-DD, N fields, M H terms, vX.Y.Z)`. This tracks regressions and sets expectations. Use `--timeout 0` for theories that exceed the default 600s.
- **Convergence testing**: After solver changes, verify error decreases at expected rate with resolution (4x for 2nd-order FD, 16x for 4th-order, machine-precision for spectral)
- **Regression detection**: Map changed files to relevant physics tests (see `/validate-physics` skill). Run those tests, not the full suite, for fast feedback.
- **t_end independence test for conversion amplification**: After measuring P_torsion/P_GR, ALWAYS verify at two different t_end values (e.g., t and 2t). If A(2t)/A(t) ≈ 1 → genuine amplification. If A(2t)/A(t) >> 1 → tachyonic instability artifact (see #238). B₀ scaling does NOT distinguish amplification from instability (growth rate is B₀-independent). IC amplitude must satisfy **h ≪ 1** (metric linearization; errors are O(h²), so h = 0.1 gives ~1% second-order corrections). Do NOT compare h directly to B₀ — they play different dimensional roles and the null result is IC-amplitude-independent since P = sin²(κB₀t/2).
- **Perturbative P regime**: For valid amplification factor A = P/P_GR, ensure P_GR ≪ 1 (equivalently κB₀D/2 ≪ 1). When P is not small, A becomes B₀-dependent and the linearized conversion measurement breaks down. Always verify P_max < 0.1 across the entire sweep before trusting A values. Choose B₀ and t_end accordingly.

## Common Pitfalls

- **Underscore constants**: `B0_peak` → `Pattern[B0, Blank[peak]]` in Mathematica. Use `Bpeak`.
- **Negative CLI values**: use `=` syntax: `--bounds="-100:100"` (not `--bounds "-100:100"`)
- **Memory size**: MEMORY.md must stay under 200 lines (excess silently truncated)
- **Implicit-dynamical sector, closure restriction, gauge certificate (GH #457/#468)**: order-0 rows can be secretly dynamical — `EquationSystem.second_order_sector` is the ONE definition (call them "implicit-dynamical sector fields"; "promoted" is reserved for ε-promotion of new DOF, #321). On the localized implicit-dynamical class the full pencil refuses at float64 (#474); `tidal simulate` then evolves the exactly closed sector excited by the IC (`restricted_spec.json` becomes the run's `spec_path`; omitted fields are ABSENT, never zero; `--no-closure-restriction` attempts the full system). Every `tidal measure` result carries a `gauge_certificate` — `flagged` means the number depends on the pinned (min-norm-zero) gauge choice; never present a flagged value as gauge-invariant. Theory and decision record: `docs/tex/pencil_engine.tex`.
- **Wolfram Exp overflow**: serializes `Exp[-x²]` as `1/E^(x²)` → Python overflow. Use `_invert_exp_denominator()`.
- **Plane-wave IC snap on periodic grids**: `--ic plane-wave --ic-wavevector k0` auto-snaps `k0` to the nearest discrete Fourier mode `2π·n/L` (clamped below Nyquist) to eliminate spectral leakage. A `Note:` is printed when the snap is significant. Off-grid `k` causes cos(k·x) to leak amplitude onto every discrete mode; for theories with tachyonic eigenvalues at some k-modes (e.g. PGT torsion with cross-coupling) this triggers spurious `SimulationDivergedError`. Sharp 0.1%-scale stability boundaries in parameter space are a signature of this — if seen, check whether the IC wavevector is off-grid. Override with `--ic-no-snap`. See `docs/tex/plane_wave_ic.tex`.
- **FV ↔ TorsionCDT dark-photon equivalence map (current convention, post-2026-04-24)**: The CDT Lagrangian in `dark_photon_plasma/theory.toml` is `L ⊃ -alpha3·I3` (note the leading minus, standard Proca convention after the 2026-04-24 sign-flip), giving spatial EOM `m² = +2·alpha3`. FV writes `L ⊃ -½·mT²·t·t` → `m² = +mT²`. Equivalence map: **`mT2 = 2·alpha3`** (same-sign). `alpha3 > 0` ↔ `mT2 > 0` = stable Proca. Confirmed bit-exact (Δ/P_max ≈ 1e-14) at campaign parameters; see `examples/dark_photon_plasma/theory.toml:18-35` and `fv_cdt_equivalence_verified.md`. **Old convention (pre-2026-04-24)**: Lagrangian was `+alpha3·I3`, giving `m² = -2·alpha3` and equivalence `mT2 = -2·alpha3` (opposite-sign). All old-convention HPC runs (28216072, 28226826, 28366464) are archived in CAMPAIGN.md as "wrong regime for the Proca dark-photon analogy" — do not reuse those results under the new convention without re-interpretation.

## Reading Equation Specifications

**Never read a spec JSON directly, and never write a fresh scan over `examples/data/`.**
Both are how six confidently-wrong diagnoses were produced (GH #401). A torsion spec is
~96,000 tokens of JSON; the commands below answer a question in a few hundred, and the
answer comes from a vetted accessor rather than from inference by eye.

| question | command | ~tokens |
| --- | --- | --- |
| what is this coefficient, and where else is it recorded? | `tidal inspect SPEC --coefficient 'h_5:identity(h_5)'` | 200 |
| which components belong together, which are temporal? | `tidal inspect SPEC --families` | 130 |
| what does one equation say? | `tidal inspect SPEC --equation h_5` (accepts `a,b` or `all`) | 760 |
| did re-derivation change the physics? | `tidal inspect OLD --diff NEW` (exit 1 = real change) | varies |

**`--diff` is for legacy re-derivations only.** It is *not* the port gate for `tidalcosmo`:
#513 has the new package emit CAMB/PSALTer conventions natively, so a byte diff would report
every intended change as a failure. The frozen legacy spec is the oracle *for the physics* --
same equations, coefficients and signs under the new naming, established by comparison and
recorded as a written mapping (#525).
| just the proven signs, whole spec | `tidal inspect SPEC --detail summary` | 4,000 |

**Corpus-level questions are already answered.** `tests/data/spec_semantics.txt` is a
committed report of families, index structure and proven sign conflicts for every spec.
Read it instead of scanning. Regenerate with `python -m scripts.spec_semantics_report`.

**Prefer the text output to `--json`** — measured, `--json` costs ~3x more and buys
nothing when you are reading rather than parsing. `--json` is for scripts.

Two properties worth relying on: a sign verdict says `unknown` rather than guessing when
it cannot be proven, and every verdict names the tactic that decided it. Only ~8% of
coefficients have a provable sign — most are free sweep parameters — so use
`--assume-positive`/`--assume-nonzero` when you have physical grounds. Python API:
`tidal.symbolic.spec_query` (accessors) and `tidal.symbolic.sign_algebra` (the sign
decisions). See `/spec` for worked recipes.

## Cosmology Program (ACTIVE — the current direction)

The project has pivoted (2026-08-29) to a **Cobaya extension**: evolve a new sector's linear
perturbations as **spectators** on a CAMB LCDM background and turn them into CMB observables
and real likelihoods. Umbrella **#488**; the operational record is `docs/COSMOLOGY_PROGRAM.md`
(read that first — it carries the decisions register, the observable ladder and the
workstreams). Eleven design documents live in `docs/cosmology/`.

- **Two packages coexist.** `tidal/` is **legacy**; `tidalcosmo/` is the new package, written
  clean beside it (strangler fig). `tidalcosmo` is a placeholder name — it renames to `tidal`
  once legacy is deleted (M7).
- **New code NEVER imports legacy `tidal/`.** No adapters, no shims. Capabilities are *ported*
  — redesigned and moved with docstrings and issue references — never wrapped. Enforced by
  `tests_cosmo/test_package_boundary.py` and `tests/test_repo_hygiene.py`.
- **Two test suites.** `tests/` (legacy) stays untouched and green; new tests go in
  `tests_cosmo/`. Both are in `testpaths`, coverage, pyright and CI.
- **Legacy stays runnable** until every capability we want is ported and verified against the
  frozen oracle (#525), then it is deleted per capability.
- **#513 changes derive conventions**: the new package emits CAMB and PSALTer conventions
  *natively* rather than translating. Gauge is an explicit named input and travels in the spec
  as metadata.
- **Wolfram is derivation-time only.** Nothing symbolic at sampling time — that is what the
  two-stage spectrum architecture (#495) and the eikonal reduction (#504) both exist to ensure.
- **D4: no HPC without explicit permission. Local only.** See `docs/hpc_workflow.md`.

Legacy documentation under `docs/tex/` correctly describes code that still exists and still
works; it is not stale and is not to be pre-emptively rewritten. It becomes wrong when the
corresponding capability is actually replaced.

## Claude Code Skills

Custom commands in `.claude/skills/` (main conversation only, not available to subagents):
- `/test [args]` — Smart-scope pytest (auto-detects relevant tests from git diff)
- `/derive <toml>` — Safe Wolfram derivation (blocks parallel runs, validates, smoke tests)
- `/validate` — Full pipeline validation with auto-fix (lint → types → spell → tests → simulate)
- `/backup` — Memory backup and MEMORY.md health check
- `/commit [message]` — Conventional commit with mandatory pre-commit testing
- `/validate-physics` — Physics regression detection (maps changed solver/measurement files to relevant tests)
- `/bump [patch|minor]` — Version bump with commit analysis (suggests level, dry-run preview, git tag)
- `/sync-docs` — Review and update all documentation for accuracy (stats, phase status, resolved issues)

## Local Literature

`literature/` contains arXiv TeX sources for frequently-cited papers (Gertsenshtein, torsion, axion-photon mixing, and the cosmology-program set). **Always check `literature/` before searching online.** Read the TeX source directly — it's faster and more reliable than web fetches. For new frequently-cited papers, download TeX via arXiv and add to `literature/<arxiv-id>/`. `literature/README.md` is generated and tracked; the curated index with per-paper notes is `docs/references.md`.

## Architecture Reference

See `docs/tex/architecture.tex` for the complete architecture reference covering: solver backends, E-L velocity form, mass/coupling matrices, Christoffel computation, background fields, gauge fixing, xAct patterns, operators, examples, and known issues. `docs/README.md` indexes every document.

See also: `docs/tex/troubleshooting.tex`, `docs/tex/background_fields.tex`, `docs/tex/constraint_fields.tex`, `docs/tex/solver_migration.tex`, `docs/tex/gauge_fixing.tex`, `docs/tex/adaptive_timestepping.tex`, `docs/tex/architecture.tex`.

## Memory Backup

Claude auto-memory files, plans, and project settings are backed up to `.claude-memory-backup/`, `.claude-plans-backup/`, and `.claude-project-backup/` (all git-ignored). On container rebuild, all are auto-restored from backup if the volume is empty. Manual sync: `bash .devcontainer/scripts/sync-claude-memory.sh backup|restore|status`.

## HPC Workflow (CSD3) — archived

**Inactive under D4** (no HPC without explicit permission; local only — the masters project
has ended, so large HPC computation now needs permission first). The full workflow is
preserved verbatim at `docs/hpc_workflow.md` and applies again once HPC work is authorized.

## Session Persistence Workaround

The VS Code Claude Code extension has a known bug where past conversations disappear from the dropdown on window reload (upstream: https://github.com/anthropics/claude-code/issues/18619). Session `.jsonl` files persist on disk but the index files are never written. To rebuild the index and restore sessions in the dropdown: `bash .devcontainer/scripts/reindex-claude-sessions.sh`. This runs automatically on container rebuild via `postCreateCommand`.
