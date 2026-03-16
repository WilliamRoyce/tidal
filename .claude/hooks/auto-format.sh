#!/bin/bash
# PostToolUse hook: auto-format Python files after Edit/Write
# Errors are NOT suppressed — ruff failures will be visible
INPUT=$(cat)
FILE=$(echo "$INPUT" | jq -r '.tool_input.file_path // .tool_input.path // empty')
[[ "$FILE" == *.py ]] || exit 0
cd /workspaces/torsion-gertsenshtein
uv run ruff format --quiet "$FILE"
uv run ruff check --fix --quiet "$FILE"
exit 0
