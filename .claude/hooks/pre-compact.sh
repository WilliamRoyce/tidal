#!/bin/bash
# PreCompact hook: inject branch/commit context before context compaction
# After compaction, Claude loses conversation history — this breadcrumb helps
cd /workspaces/torsion-gertsenshtein 2>/dev/null || exit 0
BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
RECENT=$(git log --oneline -5 2>/dev/null || echo "no commits")
MODIFIED=$(git diff --name-only HEAD 2>/dev/null | head -10)
STAGED=$(git diff --name-only --staged 2>/dev/null | head -10)

# Escape newlines for JSON
RECENT_ESC=$(echo "$RECENT" | tr '\n' '|' | sed 's/|/\\n/g')
MODIFIED_ESC=$(echo "$MODIFIED" | tr '\n' '|' | sed 's/|/\\n/g')
STAGED_ESC=$(echo "$STAGED" | tr '\n' '|' | sed 's/|/\\n/g')

cat <<EOF
{
  "additionalContext": "CONTEXT BEFORE COMPACTION:\\nBranch: $BRANCH\\nRecent commits:\\n$RECENT_ESC\\nModified files:\\n$MODIFIED_ESC\\nStaged files:\\n$STAGED_ESC\\nRemember: check CLAUDE.md and memory files for project rules after compaction."
}
EOF