---
name: sync-docs
description: Review and update all documentation for accuracy. Use when docs have drifted out of sync with the codebase, or before a release.
---

# Documentation Sync

## Current project stats
!`echo "Tests: $(uv run pytest tests/ --collect-only -q 2>&1 | tail -1)" && echo "Examples: $(ls examples/*/theory.toml 2>/dev/null | wc -l) theory.toml files" && echo "JSON specs: $(ls examples/data/*.json 2>/dev/null | wc -l)" && echo "Skills: $(ls .claude/skills/*/SKILL.md 2>/dev/null | wc -l)" && echo "Version: $(python3 -c "import re; m=re.search(r'version = \"([^\"]+)\"', open('pyproject.toml').read()); print(m.group(1) if m else '?')")"

## Open ROADMAP items
!`grep -n "🔄\|TODO\|BLOCKED\|In Progress" docs/ROADMAP.md 2>/dev/null | head -15`

## Recent commits (for context on what changed)
!`git log --oneline -15`

## Instructions

### Step 1 — Identify stale documentation
Read through key docs and compare against current codebase state:
- `docs/ROADMAP.md` — Are phase statuses current? Any resolved issues still marked open?
- `docs/NEXT_PHASES.md` — Are completed phases marked done? Are stats accurate?
- `README.md` — Are test counts, example counts, feature lists accurate?
- `CONTRIBUTING.md` — Are test counts accurate?
- `docs/tex/troubleshooting.tex` — Are solved issues still listed as unresolved?
- Active feature checklists — Are completed items marked done?

### Step 2 — Fix discrepancies
Update each file with current information. For numeric stats (test counts, example counts), use the actual numbers from the stats section above.

For substantive content:
- Mark resolved ROADMAP issues with how they were resolved
- Update phase completion status
- Remove or mark solved troubleshooting entries
- Update benchmark numbers if they've changed

### Step 3 — Commit
Stage and commit all documentation updates:
```bash
git add README.md CONTRIBUTING.md docs/ && git commit -m "docs: sync documentation with current project state"
```

### Step 4 — Update MEMORY.md
If test count, example count, or major feature status changed, update the corresponding lines in MEMORY.md.
