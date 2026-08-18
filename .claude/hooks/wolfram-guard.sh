#!/bin/bash
# PreToolUse hook: block parallel wolframscript (single Wolfram Engine license)
# Deny mechanism: exit code 2 = block the action
INPUT=$(cat)
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty')

# Only check commands that actually invoke wolframscript
echo "$COMMAND" | grep -qE '(wolframscript|uv run tidal derive)' || exit 0

# --- False-positive exclusions (GH #400) ---------------------------------
# The match above is on command TEXT, so it also catches commands that merely
# mention wolframscript, or tidal derive invocations that never start it.
# Those consume no license and must not be blocked.

# 1. Read-only inspections that happen to contain the pattern, e.g.
#    `pgrep -f wolframscript` — which is exactly what you want to run WHILE a
#    derivation is in flight.  Match on the first word of the command.
# head -1 first: $COMMAND may be multi-line (heredocs, continuations), and cut
# would otherwise emit one field per line rather than the leading word.
FIRST_WORD=$(printf '%s' "$COMMAND" | head -1 | sed -E 's/^[[:space:]]*//' | cut -d' ' -f1)
FIRST_WORD=$(basename "$FIRST_WORD" 2>/dev/null || printf '%s' "$FIRST_WORD")
case "$FIRST_WORD" in
  # read-only inspection
  pgrep|pkill|ps|grep|rg|egrep|fgrep|ag|awk|sed|cat|head|tail|wc|less|ls|find|\
  diff|cmp|jq|sort|uniq|tee|file|stat)
    exit 0
    ;;
  # tools that never invoke Wolfram but routinely quote its name — writing an
  # issue comment or commit message about wolframscript must not be blocked
  gh|git|echo|printf|cp|mv|mkdir|touch)
    exit 0
    ;;
esac

# 2. `tidal derive --dry-run` only prints the generated script and returns
#    before any wolframscript call; --help likewise.  Neither takes the license.
echo "$COMMAND" | grep -qE '(^|[[:space:]])--(dry-run|help)([[:space:]]|$)' && exit 0

# Block if wolframscript is already running.
#
# Use -x (match the process NAME) rather than -f (match the whole command
# line).  With -f, any shell whose command line merely *mentions*
# wolframscript matches — including the very shell running the command being
# checked — so the guard fired spuriously on commands that were only talking
# about wolframscript.  See GH #400.
if pgrep -x wolframscript > /dev/null 2>&1 || pgrep -x WolframKernel > /dev/null 2>&1; then
  echo "BLOCKED: wolframscript already running (single Wolfram Engine license). Wait for it to finish or kill it first." >&2
  exit 2
fi
exit 0
