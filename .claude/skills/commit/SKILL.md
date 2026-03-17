---
name: commit
description: Conventional commit with mandatory pre-commit testing and auto-formatting. Use when ready to commit changes.
---

# Conventional Commit

## Git status
!`cd /workspaces/torsion-gertsenshtein && git status --short`

## Changed files
!`cd /workspaces/torsion-gertsenshtein && (git diff --name-only HEAD 2>/dev/null; git diff --name-only --staged 2>/dev/null; git ls-files --others --exclude-standard 2>/dev/null) | sort -u`

## Recent commits (for style matching)
!`cd /workspaces/torsion-gertsenshtein && git log --oneline -5`

## Version status
!`cd /workspaces/torsion-gertsenshtein && echo "v$(python3 -c "import re; print(re.search(r'^version\s*=\"([^\"]+)\"', open('pyproject.toml').read(), re.MULTILINE).group(1))" 2>/dev/null) — $(git log --oneline $(git log --all --grep='chore: bump version' --format='%H' -1 2>/dev/null || git rev-list --max-parents=0 HEAD)..HEAD 2>/dev/null | wc -l) commits since last bump"`

## Instructions

### Step 1 — Run relevant tests
Map changed source files to test files:
- `tidal/solver/X.py` → `tests/test_solver_X.py`
- `tidal/cli/_X.py` → `tests/test_cli.py`, `tests/test_cli_parsing.py`
- `tidal/measurement/` → `tests/test_measurement.py`, `tests/test_new_measurements.py`
- `tidal/symbolic/` → `tests/test_json_loader.py`
- If unsure → `uv run pytest tests/ -x -q`

If ANY test fails: **STOP. Do not commit.** Report the failure and suggest fixes.

### Step 2 — Auto-format all changed files
```bash
uv run ruff check --fix
uv run ruff format
```
Stage any files that were auto-fixed.
Note: The PostToolUse hook auto-formats files Claude edits, but this step catches files the user edited manually outside Claude.

### Step 3 — Draft commit message
Use conventional commit format:
- `feat:` — new feature
- `fix:` — bug fix
- `refactor:` — restructuring without behavior change
- `test:` — test changes only
- `docs:` — documentation only
- `perf:` — performance improvement

First line under 72 characters. Add a body paragraph for non-trivial changes explaining "why".
Include references when relevant (paper citations, issue numbers, algorithm names).

### Step 4 — Commit
Stage SPECIFIC files by name (never `git add -A` or `git add .`).
Separate unrelated changes into distinct commits.

### Step 5 — Version bump (if task is complete)
After committing, evaluate whether to bump the version:
- **Bump patch** (0.0.x): if committing a `fix:`, `refactor:`, `perf:`, `docs:` change that completes a task
- **Bump minor** (0.x.0): if committing a `feat:` that completes a new feature
- **Skip bump if**: this is a WIP commit, part of a multi-commit feature still in progress, or tests are not fully passing

If bumping:
```bash
python scripts/bump_version.py --{level} --commit --allow-dirty
```
NEVER bump the major version automatically — major bumps require explicit user request.

### CRITICAL RULES
- **NO Co-Authored-By trailer** — never add it
- **NO committing .env, credentials, or large binaries**
- If `$ARGUMENTS` is provided, use it as a hint for the commit message
- If changes span unrelated areas, make separate commits for each
