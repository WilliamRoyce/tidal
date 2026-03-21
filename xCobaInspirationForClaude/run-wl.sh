#!/bin/sh
# Run a Wolfram Language script with xperm cleanup.
# Usage: ./run-wl.sh Script.m [logfile]

if [ -z "$1" ]; then
    echo "Usage: $0 <script.m> [logfile]" >&2
    exit 1
fi

SCRIPT="$1"
LOGFILE="${2:-/tmp/$(basename "$SCRIPT" .m).log}"

wolfram -script "$SCRIPT" > "$LOGFILE" 2>&1
pgrep -x xperm.linux.64-bit | xargs -r kill 2>/dev/null
true
