---
name: validate
description: Full pipeline validation with auto-fix for lint, type, and spell errors. Use before merging or when the user says "validate".
---

# Full Pipeline Validation (Auto-Fix)

## Current status
!`cd /workspaces/torsion-gertsenshtein && echo "=== Lint ===" && uv run ruff check --output-format concise 2>&1 | tail -10 && echo "=== Type ===" && uv run pyright --outputjson 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(f'{len(d.get(\"generalDiagnostics\",[]))} type errors')" 2>&1 || echo "pyright needed"`

## Instructions

Run each step sequentially. Stop and fix issues before proceeding to the next step.

### Step 1 — Lint (auto-fix)
```bash
uv run ruff check --fix
uv run ruff format
```
Report what was auto-fixed. If any unfixable errors remain, read the files and fix them manually.

### Step 2 — Type errors (auto-fix)
Run `uv run pyright`. For EACH error reported:
- Read the offending file and understand the context
- Fix the type issue (add annotations, handle Optional, fix imports, correct return types)
- Only ask the user if the fix is genuinely ambiguous

### Step 3 — Spelling (auto-fix)
Run `uv run cspell "tidal/**/*.py" "tests/**/*.py"`. For each flagged word:
- Physics/math domain terms (xAct, Christoffel, Lagrangian, Hamiltonian, Fornberg, Yoshida, Boccaletti, eigendecomp, etc.) → add to `.cspell.json` words list
- Genuine typos → fix in the code

### Step 4 — Tests
```bash
uv run pytest tests/ -x -q
```
If any test fails: analyze the root cause, fix it, and re-run. Continue until all tests pass.
The test suite includes 80+ simulation integration tests with conservation checks at 1e-6 threshold — no separate smoke simulation needed.

### Step 5 — Documentation accuracy
Check key documentation for obvious staleness:
```bash
echo "Actual tests: $(uv run pytest tests/ --collect-only -q 2>&1 | tail -1)"
echo "Actual examples: $(ls examples/*/theory.toml 2>/dev/null | wc -l)"
grep -n "Python tests" README.md docs/ROADMAP.md docs/NEXT_PHASES.md 2>/dev/null
```
If discrepancies found, warn and suggest running `/sync-docs`.

### Step 6 — Report
Summary table of everything found and fixed:
- Lint: N issues auto-fixed, M remaining
- Types: N errors fixed
- Spelling: N words added to dictionary, M typos fixed
- Tests: N passed, M failed (with details)
- Docs: in sync / N discrepancies found

### Step 7 — File issues for problems found
If validation revealed notable problems (whether fixed or not), create GitHub issues to build a trail:
- Issues you fixed during validation → create and close immediately with the fix details
- Issues that couldn't be auto-fixed (architectural, design, persistent failures) → create and leave open
Check for duplicates first: `gh issue list -S "keyword"`.
