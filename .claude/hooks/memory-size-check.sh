#!/bin/bash
# PostToolUse hook: warn if MEMORY.md exceeds 200-line limit
INPUT=$(cat)
FILE=$(echo "$INPUT" | jq -r '.tool_input.file_path // .tool_input.path // empty')
[[ "$FILE" == */MEMORY.md ]] || exit 0
MEMFILE="$HOME/.claude/projects/-workspaces-torsion-gertsenshtein/memory/MEMORY.md"
[ -f "$MEMFILE" ] || exit 0
LINES=$(wc -l < "$MEMFILE" 2>/dev/null) || LINES=0
if [ "$LINES" -gt 200 ]; then
  echo "WARNING: MEMORY.md is $LINES lines (limit: 200). Content past line 200 is silently truncated. Move detail to topic files." >&2
fi
exit 0