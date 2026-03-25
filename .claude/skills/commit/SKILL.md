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

## Docs mentioning changed components
!`cd /workspaces/torsion-gertsenshtein && for f in $(git diff --name-only HEAD 2>/dev/null | head -5); do base=$(basename "$f" .py); grep -rl --include="*.tex" --include="*.md" "$base" docs/ 2>/dev/null; done | sort -u | head -10`

## Related open issues
!`cd /workspaces/torsion-gertsenshtein && gh issue list --limit 10 --json number,title --jq '.[] | "#\(.number): \(.title)"' 2>/dev/null`

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

### Step 5 — Version bump
After committing, bump the version. Default to bumping — only skip if you have more code changes to make in this same sitting:
- **Bump patch** (0.0.x): for `fix:`, `refactor:`, `perf:`, `docs:`, `test:`, `revert:`, `style:`, `chore:` commits
- **Bump minor** (0.x.0): for `feat:` commits that add new functionality
- **Skip ONLY if**: you are about to make another commit immediately (same sitting, same task)

If bumping:
```bash
python scripts/bump_version.py --{level} --commit --allow-dirty
```
NEVER bump the major version automatically — major bumps require explicit user request.

### Step 6 — Update documentation (if task is complete)
After committing a completed feature or fix, check if documentation needs updating.

**How to find relevant docs:** Search `docs/` for mentions of the component you changed:
```bash
grep -rl "keyword" docs/
```
Also check the "Docs mentioning changed components" section above.

**Documentation structure:**
- Technical docs live in `docs/tex/*.tex` — these are the primary documentation
- Project management (roadmaps, checklists) lives in `docs/*.md`
- See `docs/README.md` for the full documentation index

**Common patterns:**
- Resolved an issue? → Update status in `docs/ROADMAP.md` or `docs/NEXT_PHASES.md`
- Completed a checklist item? → Mark done in the relevant `docs/*checklist*.md`
- Changed performance? → Update benchmark tables in whichever `docs/tex/*.tex` covers that subsystem
- Discovered new error pattern? → Add to `docs/tex/troubleshooting.tex`
- Changed algorithm/architecture? → Update the `docs/tex/*.tex` that describes that component
- Run `/sync-docs` after major features to check for drift

If docs need updating, make the changes and commit separately:
```bash
git add docs/ && git commit -m "docs: update {topic} after {feature/fix}"
```

**Skip if**: this is a WIP commit, the change is trivial, or no docs mention the changed component.

### Step 7 — Create issues for anything notable discovered
During this work, if you encountered any of the following, create a GitHub issue to build a searchable trail:
- Bugs or unexpected behavior (even if you fixed them)
- Performance issues worth investigating
- Missing features that would have helped
- Technical debt or inconsistencies
- Design decisions worth recording
- Test gaps or coverage holes

Check for duplicates first: `gh issue list -S "keyword"`.

For things you **fixed in this commit**:
```bash
gh issue create --title "Short title" --body "..." --label "bug"
gh issue close N -c "Fixed in <commit-hash>"
```

For things **still open**:
```bash
gh issue create --title "Short title" --body "..." --label "enhancement"
```

**Skip if**: truly trivial (typo/formatting) or a duplicate already exists.
**NEVER include any "Generated with Claude Code" footer or attribution in issue bodies.**

### CRITICAL RULES
- **NO Co-Authored-By trailer** — never add it
- **NO committing .env, credentials, or large binaries**
- **Only stage files YOU changed** — before `git add`, verify each file's diff is from YOUR work in this session. Never blindly stage files that may have been modified by parallel agents or worktrees. Check `git diff <file>` if unsure.
- If `$ARGUMENTS` is provided, use it as a hint for the commit message
- If changes span unrelated areas, make separate commits for each
