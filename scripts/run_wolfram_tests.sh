#!/bin/bash
# Run all Wolfram unit tests
#
# Usage: ./scripts/run_wolfram_tests.sh
#
# Runs each test file in tests/wolfram/ and reports summary.
# Exit code: 0 if all pass, 1 if any fail.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Running Wolfram Tests ==="
echo ""

# Discover all test files
TESTS=(
    "tests/wolfram/test_euler_lagrange.wls"
    "tests/wolfram/test_common_utilities.wls"
    "tests/wolfram/test_export_json.wls"
    "tests/wolfram/test_component_decompose.wls"
)

PASSED=0
FAILED=0
FAILED_TESTS=()

for test in "${TESTS[@]}"; do
    TEST_PATH="$PROJECT_ROOT/$test"

    if [ ! -f "$TEST_PATH" ]; then
        echo "SKIP: $test (file not found)"
        continue
    fi

    echo "Running: $test"
    if wolframscript -file "$TEST_PATH"; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
        FAILED_TESTS+=("$test")
    fi
    echo ""
done

echo "=== Wolfram Test Summary ==="
echo "Passed: $PASSED"
echo "Failed: $FAILED"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Failed tests:"
    for t in "${FAILED_TESTS[@]}"; do
        echo "  - $t"
    done
    exit 1
fi

echo ""
echo "*** ALL WOLFRAM TESTS PASSED ***"
exit 0
