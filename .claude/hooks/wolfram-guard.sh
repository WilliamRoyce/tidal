#!/bin/bash
# PreToolUse hook: block parallel wolframscript (single Wolfram Engine license)
# Deny mechanism: exit code 2 = block the action
INPUT=$(cat)
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty')

# Only check commands that actually invoke wolframscript
echo "$COMMAND" | grep -qE '(wolframscript|uv run tidal derive)' || exit 0

# Block if wolframscript is already running
if pgrep -f wolframscript > /dev/null 2>&1; then
  echo "BLOCKED: wolframscript already running (single Wolfram Engine license). Wait for it to finish or kill it first." >&2
  exit 2
fi
exit 0
