---
name: validate-physics
description: Detect physics regressions from recent code changes. Maps changed solver/measurement files to relevant physics tests and runs them. Use after modifying solver backends, measurement code, or operator implementations.
---

# Physics Regression Detection

## Recent changes to physics code
!`cd /workspaces/torsion-gertsenshtein && git diff --name-only HEAD~5 -- tidal/solver/ tidal/measurement/ tidal/symbolic/ 2>/dev/null | sort -u`

## Instructions

### Step 1 — Map changes to relevant tests
Identify which physics tests cover the changed code:
- `tidal/solver/modal.py` → `tests/test_solver_modal.py`
- `tidal/solver/ida.py` or `cvode.py` → `tests/test_solver_ida.py`, `tests/test_solver_cvode.py`
- `tidal/solver/leapfrog.py` → `tests/test_solver_leapfrog.py`
- `tidal/solver/operators.py` → `tests/test_operators.py`
- `tidal/solver/fields.py` or `state.py` → `tests/test_solver_fields.py`
- `tidal/solver/rhs.py` or `coefficients.py` → `tests/test_solver_rhs.py`
- `tidal/measurement/` → `tests/test_measurement.py`, `tests/test_new_measurements.py`
- `tidal/symbolic/` → `tests/test_json_loader.py`
- If changes span many files → run full suite

### Step 2 — Run the mapped tests
```bash
uv run pytest <mapped_test_files> -x -v
```
Use `-v` (verbose) instead of `-q` — physics regressions need detailed output to diagnose.

### Step 3 — Analyze any failures
If tests fail:
- Read the test to understand what physics property it checks
- Read the changed code to understand what was modified
- Determine if the failure is a genuine regression or an expected change
- Fix regressions; update test expectations only if the physics changed intentionally

### Step 4 — Report
- Which files changed and which tests were run
- Pass/fail results with failure details
- Whether failures are regressions vs expected changes

### Step 5 — File issues for regressions
If physics tests revealed regressions or notable behavior:
- Regressions you fixed → create issue with label `validation`, then close with fix commit
- Unclear root cause or needs investigation → create issue with labels `validation` + `needs-investigation`, leave open
- Include: affected example, observed vs expected behavior, relevant file paths
Check for duplicates first: `gh issue list -S "keyword"`.