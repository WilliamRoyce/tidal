#!/bin/bash
# Check Wolfram modules for basic syntax errors
#
# Usage: ./scripts/lint_wolfram.sh
#
# Loads each module to verify it has no syntax errors.
# Does not run any tests, just validates modules can be loaded.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Wolfram Module Syntax Check ==="
echo ""

MODULES=(
    "tidal/wolfram/CommonUtilities.wl"
    "tidal/wolfram/ExportJSON.wl"
    "tidal/wolfram/ComponentDecompose.wl"
    "tidal/wolfram/EulerLagrange.wl"
    "tidal/wolfram/Linearize.wl"
    "tidal/wolfram/LagrangianPipeline.wl"
)

PASSED=0
FAILED=0
FAILED_MODULES=()

for module in "${MODULES[@]}"; do
    MODULE_PATH="$PROJECT_ROOT/$module"

    if [ ! -f "$MODULE_PATH" ]; then
        echo "SKIP: $module (file not found)"
        continue
    fi

    echo -n "Checking: $module ... "
    if wolframscript -c "Quiet[Get[\"$MODULE_PATH\"]]; Print[\"OK\"]" 2>/dev/null | grep -q "OK"; then
        echo "OK"
        ((PASSED++))
    else
        echo "FAILED"
        ((FAILED++))
        FAILED_MODULES+=("$module")
    fi
done

echo ""
echo "=== Syntax Check Summary ==="
echo "Passed: $PASSED"
echo "Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Failed modules:"
    for m in "${FAILED_MODULES[@]}"; do
        echo "  - $m"
    done
    exit 1
fi

echo ""
echo "*** All modules loaded successfully ***"
exit 0
