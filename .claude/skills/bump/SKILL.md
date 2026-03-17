---
name: bump
description: Version bump with conventional-commit analysis. Suggests patch/minor based on commit history since last bump. Use when ready to release a version.
---

# Version Bump

## Current version
!`python3 -c "import re; print(re.search(r'^version\s*=\"([^\"]+)\"', open('pyproject.toml').read(), re.MULTILINE).group(1))"`

## Commits since last version bump
!`git log --oneline $(git log --all --grep='chore: bump version' --format='%H' -1 2>/dev/null || git rev-list --max-parents=0 HEAD)..HEAD 2>/dev/null | head -40`

## Commit type breakdown
!`git log --format="%s" $(git log --all --grep='chore: bump version' --format='%H' -1 2>/dev/null || git rev-list --max-parents=0 HEAD)..HEAD 2>/dev/null | grep -oP '^\w+(?=:)' | sort | uniq -c | sort -rn`

## Instructions

### Step 1 — Determine bump level
If `$ARGUMENTS` contains `patch` or `minor`, use that directly.

Otherwise, analyze the commit type breakdown above:
- Any `feat:` commits → **minor** (0.x.0)
- Only `fix:`/`refactor:`/`docs:`/`test:`/`perf:`/`style:`/`chore:` → **patch** (0.0.x)

**NEVER bump the major version.** Only `--patch` and `--minor` are available.
Major bumps require explicit version: `python scripts/bump_version.py X.0.0 --commit`.

### Step 1.5 — Sync documentation
Before bumping, ensure docs are current. Check key stats (test counts, example counts) and ROADMAP status. If discrepancies found, fix them so they're included in the bump.

### Step 2 — Dry-run preview
```bash
python scripts/bump_version.py --{level} --dry-run
```
Show: current → new version, commit count, type breakdown, files affected.

### Step 3 — Execute bump
```bash
python scripts/bump_version.py --{level} --commit --allow-dirty
```

### Step 4 — Report
- Old → new version, git tag created
- Commit count and type breakdown
- Remind: `git push && git push --tags` when ready