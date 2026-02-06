#!/bin/bash
# Full test suite: Python tests + Wolfram tests
#
# Usage: ./scripts/full_test.sh
#
# Runs both Python (pytest) and Wolfram test suites.
# Good to run before committing changes.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

echo "=== Full Test Suite ==="
echo ""

# Track overall status
PYTHON_OK=0
WOLFRAM_OK=0

# Python tests
echo ">>> Python Tests (pytest)"
echo ""
if uv run pytest tests/ -v; then
    PYTHON_OK=1
    echo ""
    echo "Python tests: PASSED"
else
    echo ""
    echo "Python tests: FAILED"
fi
echo ""

# Wolfram tests
echo ">>> Wolfram Tests"
echo ""
if "$SCRIPT_DIR/run_wolfram_tests.sh"; then
    WOLFRAM_OK=1
    echo ""
    echo "Wolfram tests: PASSED"
else
    echo ""
    echo "Wolfram tests: FAILED"
fi
echo ""

# Summary
echo "=== Full Test Summary ==="
echo "Python:  $([ $PYTHON_OK -eq 1 ] && echo 'PASSED' || echo 'FAILED')"
echo "Wolfram: $([ $WOLFRAM_OK -eq 1 ] && echo 'PASSED' || echo 'FAILED')"
echo ""

if [ $PYTHON_OK -eq 1 ] && [ $WOLFRAM_OK -eq 1 ]; then
    echo "*** ALL TESTS PASSED ***"
    exit 0
else
    echo "*** SOME TESTS FAILED ***"
    exit 1
fi
