#!/bin/bash
# PostToolUse hook: warn if MEMORY.md exceeds 200-line limit
INPUT=$(cat)
FILE=$(echo "$INPUT" | jq -r '.tool_input.file_path // .tool_input.path // empty')
[[ "$FILE" == */MEMORY.md ]] || exit 0
# Derive the Claude project slug: the absolute project path with every
# non-alphanumeric character replaced by "-".
PROJECT_ROOT="$(cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}" && pwd)"
PROJECT_SLUG="$(printf '%s' "$PROJECT_ROOT" | sed 's/[^a-zA-Z0-9]/-/g')"
MEMFILE="$HOME/.claude/projects/$PROJECT_SLUG/memory/MEMORY.md"
[ -f "$MEMFILE" ] || exit 0
LINES=$(wc -l < "$MEMFILE" 2>/dev/null) || LINES=0
if [ "$LINES" -gt 200 ]; then
  echo "WARNING: MEMORY.md is $LINES lines (limit: 200). Content past line 200 is silently truncated. Move detail to topic files." >&2
fi
exit 0