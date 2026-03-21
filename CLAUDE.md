# TIDAL: Tensor Integration and Derivation for Any Lagrangian

Symbolic physics pipeline: Lagrangian (xAct/Mathematica) -> JSON -> native PDE solver (SUNDIALS IDA/CVODE, leapfrog, scipy) with numpy spatial operators. All PDEs derive from the Lagrangian via symbolic computation -- never manually hardcode equations. Operates exclusively in the linearized regime (quadratic Lagrangians, linear PDEs).

## Project Structure

- `tidal/wolfram/` -- Wolfram pipeline modules (EulerLagrange.wl, ComponentDecompose.wl, ExportJSON.wl, CommonUtilities.wl, GaugeFix.wl)
- `tidal/solver/` -- PDE solver backends (ida.py, cvode.py, leapfrog.py, fields.py, operators.py, grid.py, coefficients.py, rhs.py, state.py, constraint_solve.py)
- `tidal/symbolic/` -- Python symbolic pipeline (_derive.py, json_loader.py)
- `tidal/cli/` -- CLI entry points (9 subcommands: derive, simulate, measure, inspect, list, validate, plot, sweep, analyze)
- `tidal/measurement/` -- Physics measurements (energy, conversion, mixing, spectra)
- `examples/` -- 20 physics examples (1+1D through 3+1D), each with theory.toml + .wls + data/*.json
- `tests/` -- ~1,700 Python tests + ~115 Wolfram tests
- `docs/` -- Architecture docs (MEMORY.md is the main reference)

## Key Commands

- `uv run pytest tests/ -x -q` -- Run Python tests
- `./scripts/full_test.sh` -- Full test suite (Python + Wolfram)
- `uv run tidal derive examples/<name>/theory.toml` -- Derive PDEs from Lagrangian
- `uv run tidal simulate examples/data/<name>.json` -- Run simulation
- `uv run tidal sweep examples/data/<name>.json --sweep "param=start:stop:N" --measure conversion --output sweep_out` -- Run parameter sweep
- `uv run ruff check` / `uv run ruff format` -- Lint / format
- `uv run pyright` -- Type checking

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

## Workflow Rules

- **After completing any code change**, run relevant tests before moving on. Source→test mapping: `tidal/solver/X.py` → `tests/test_solver_X.py`, `tidal/cli/_X.py` → `tests/test_cli.py`, `tidal/measurement/` → `tests/test_measurement.py`. Unsure → full suite: `uv run pytest tests/ -x -q`
- **After completing a feature/fix**, commit promptly with conventional format (feat:/fix:/refactor:/test:/docs:). No Co-Authored-By trailer. Separate unrelated changes into distinct commits.
- **Fix lint/type/spell errors immediately** — `uv run ruff check --fix && uv run ruff format` after code changes. Fix pyright errors. Add domain terms to `.cspell.json`, fix genuine typos.
- **Wolfram pipeline integrity**: ALL symbolic processing stays in Wolfram — never post-process equations in Python. Never skip/bypass the canonical pipeline; fix root causes.
- **Run long commands in background**: Use `run_in_background: true` on the Bash tool for any command that takes more than a few seconds — derivations (`tidal derive`), simulations (`tidal simulate`), sweeps (`tidal sweep`), and full test suites (`pytest tests/`). Continue other work while waiting; you'll be notified on completion. Do NOT poll or sleep.
- **Only ONE wolframscript at a time** — single engine license. NEVER run `tidal derive` in parallel.
- **Use minimal test theories** (scalar_field, coupled_scalars) before expensive derivations.
- **Negative energies** may be physical with (-,+,+,+) metric convention — don't "fix" without understanding the physics.
- **Before context compaction**, update all relevant docs and memory files.
- **Version bump after completing work**: After committing a completed feature/fix (all tests pass, no remaining tasks), bump the version: `--patch` for fixes/small changes, `--minor` for new features. Skip if mid-feature or WIP. NEVER bump the major version automatically. Use `python scripts/bump_version.py --{level} --commit --allow-dirty`. Default to bumping — only skip if you are about to make another commit immediately in the same sitting.
- **Update documentation after completing work**: After committing a feature/fix, identify and update relevant docs. To find which docs are affected, search `docs/` for mentions of the changed component/feature (e.g., `grep -rl "Jacobian\|sparse" docs/`). Common update patterns: phase/issue status changes → `docs/ROADMAP.md`, `docs/NEXT_PHASES.md`; implementation substep done → active checklist in `docs/`; performance changed → whichever `docs/tex/*.tex` has benchmark tables for that subsystem; new error pattern → `docs/tex/troubleshooting.tex`; algorithm/architecture changed → whichever `docs/tex/*.tex` describes that component. See `docs/README.md` for the documentation index. Commit doc updates separately: `docs: update {topic} documentation`.
- **Only commit YOUR changes**: Before staging files, verify each file's diff contains changes YOU made in THIS session. Never stage files modified by parallel agent sessions or other worktrees. If `git status` shows unexpected modified files, check `git diff <file>` before staging.
- **Create GitHub issues proactively**: When you encounter bugs, improvement opportunities, technical debt, or notable discoveries during work, create a GitHub issue via `gh issue create` to build a searchable trail. This applies even for things you fix immediately — create the issue, then close it with the fix commit (`gh issue close N -c "Fixed in <commit>"`) so there's a record of what was found and how it was resolved. Always check for duplicates first: `gh issue list -S "keyword"`. Use appropriate labels from the existing set (bug, enhancement, documentation, validation, etc.). Include: clear title, context, relevant file paths, and why it matters. Skip only if: truly trivial (typo, formatting) or a duplicate already exists. **NEVER include any "Generated with Claude Code" footer or attribution in issue bodies.**

## Physics Coding Patterns

- **Specify success criteria before coding**: "Modal solver must agree with CVODE to RMS < 1%" — not just "implement modal solver". Include quantitative thresholds.
- **Wolfram derivations**: Read an existing .wls template first, generate new by modifying template, review diff against template before running wolframscript
- **After derivation**: Verify JSON has `canonical.hamiltonian_terms` — without this, all energy measurements fail silently. Run `tidal validate <json> --stability`.
- **Convergence testing**: After solver changes, verify error decreases at expected rate with resolution (4x for 2nd-order FD, 16x for 4th-order, machine-precision for spectral)
- **Regression detection**: Map changed files to relevant physics tests (see `/validate-physics` skill). Run those tests, not the full suite, for fast feedback.

## Common Pitfalls

- **Underscore constants**: `B0_peak` → `Pattern[B0, Blank[peak]]` in Mathematica. Use `Bpeak`.
- **Negative CLI values**: use `=` syntax: `--bounds="-100:100"` (not `--bounds "-100:100"`)
- **Memory size**: MEMORY.md must stay under 200 lines (excess silently truncated)
- **Wolfram Exp overflow**: serializes `Exp[-x²]` as `1/E^(x²)` → Python overflow. Use `_invert_exp_denominator()`.

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

`Literature/` contains arXiv TeX sources for frequently-cited papers (Gertsenshtein, torsion, axion-photon mixing). **Always check `Literature/` before searching online.** Read the TeX source directly — it's faster and more reliable than web fetches. For new frequently-cited papers, download TeX via arXiv and add to `Literature/<arxiv-id>/`.

## Architecture Reference

See `docs/MEMORY.md` for the complete architecture reference covering: solver backends, E-L velocity form, mass/coupling matrices, Christoffel computation, background fields, gauge fixing, xAct patterns, operators, examples, and known issues.

See also: `docs/tex/troubleshooting.tex`, `docs/tex/background_fields.tex`, `docs/tex/constraint_fields.tex`, `docs/tex/solver_migration.tex`, `docs/tex/gauge_fixing.tex`, `docs/tex/adaptive_timestepping.tex`, `docs/tex/architecture.tex`.

## Memory Backup

Claude auto-memory files, plans, and project settings are backed up to `.claude-memory-backup/`, `.claude-plans-backup/`, and `.claude-project-backup/` (all git-ignored). On container rebuild, all are auto-restored from backup if the volume is empty. Manual sync: `bash .devcontainer/scripts/sync-claude-memory.sh backup|restore|status`.

## Session Persistence Workaround

The VS Code Claude Code extension has a known bug where past conversations disappear from the dropdown on window reload (upstream: https://github.com/anthropics/claude-code/issues/18619). Session `.jsonl` files persist on disk but the index files are never written. To rebuild the index and restore sessions in the dropdown: `bash .devcontainer/scripts/reindex-claude-sessions.sh`. This runs automatically on container rebuild via `postCreateCommand`.
