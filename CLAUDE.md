# TIDAL: Tensor Integration and Derivation for Any Lagrangian

Symbolic physics pipeline: Lagrangian (xAct/Mathematica) -> JSON -> native PDE solver (SUNDIALS IDA/CVODE, leapfrog, scipy) with numpy spatial operators. All PDEs derive from the Lagrangian via symbolic computation -- never manually hardcode equations. Operates exclusively in the linearized regime (quadratic Lagrangians, linear PDEs).

## Project Structure

- `tidal/wolfram/` -- Wolfram pipeline modules (EulerLagrange.wl, ComponentDecompose.wl, ExportJSON.wl, CommonUtilities.wl, GaugeFix.wl)
- `tidal/solver/` -- PDE solver backends (ida.py, cvode.py, leapfrog.py, fields.py, operators.py, grid.py, coefficients.py, rhs.py, state.py, constraint_solve.py)
- `tidal/symbolic/` -- Python symbolic pipeline (_derive.py, json_loader.py)
- `tidal/cli/` -- CLI entry points (9 subcommands: derive, simulate, measure, inspect, list, validate, plot, sweep, analyze)
- `tidal/measurement/` -- Physics measurements (energy, conversion, mixing, spectra)
- `examples/` -- 25 physics examples (1+1D through 3+1D), each with theory.toml + .wls + data/*.json
- `tests/` -- ~1,343 Python tests + ~115 Wolfram tests
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
- Cross-field decomposition requires passing `additionalFields` to `DecomposeToComponents`
- Background fields declared via `[[background_fields]]` TOML section
- Gauge fixing via `[[gauge]]` TOML section (presets: Lorenz, de Donder, Coulomb, temporal, axial)
- Velocity naming: v_{field_name} (e.g., v_phi_0, v_A_1) — E-L velocity form, not canonical momenta

## Architecture Reference

See `docs/MEMORY.md` for the complete architecture reference covering: solver backends, E-L velocity form, mass/coupling matrices, Christoffel computation, background fields, gauge fixing, xAct patterns, operators, examples, and known issues.

See also: `docs/troubleshooting.md`, `docs/background_fields.md`, `docs/constraint_fields.md`, `docs/solver_migration.md`, `docs/gauge_fixing.md`, `docs/adaptive_timestepping.md`, `docs/architecture/README.md`.

## Memory Backup

Claude auto-memory files are backed up to `.claude-memory-backup/` (git-ignored). On container rebuild, memory is auto-restored from this backup if the volume is empty. Manual sync: `bash .devcontainer/scripts/sync-claude-memory.sh backup|restore|status`.

## Session Persistence Workaround

The VS Code Claude Code extension has a known bug where past conversations disappear from the dropdown on window reload (upstream: https://github.com/anthropics/claude-code/issues/18619). Session `.jsonl` files persist on disk but the index files are never written. To rebuild the index and restore sessions in the dropdown: `bash .devcontainer/scripts/reindex-claude-sessions.sh`. This runs automatically on container rebuild via `postCreateCommand`.
