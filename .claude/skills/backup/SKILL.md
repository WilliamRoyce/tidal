---
name: backup
description: Backup Claude memory, plans, and project settings. Check MEMORY.md health. Use after significant work or before ending a session.
---

# Memory Backup & Health Check

## Current MEMORY.md size
!`wc -l "$HOME/.claude/projects/$(printf '%s' "$PWD" | sed 's/[^a-zA-Z0-9]/-/g')/memory/MEMORY.md" 2>/dev/null`

## Instructions

1. Run the backup script:
```bash
bash .devcontainer/scripts/sync-claude-memory.sh backup
```

2. Show status:
```bash
bash .devcontainer/scripts/sync-claude-memory.sh status
```

3. Check MEMORY.md line count. If over 200 lines, WARN:
   "MEMORY.md is N lines (limit: 200). Content past line 200 is silently truncated every session. Move detailed content into topic files."

4. Report: number of memory files backed up, plan files, last sync timestamp, MEMORY.md health.
