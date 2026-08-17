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

# Discover all test files.
#
# Previously a hardcoded list, which silently omitted four suites that existed
# on disk but were never executed (test_canonicalize_exp_fraction,
# test_gauge_fix, test_lps, test_metric_component_values).  Globbing means a
# new tests/wolfram/test_*.wls is picked up automatically and cannot be
# forgotten.  test_harness.wl is a library, not a suite, and is excluded by
# the test_*.wls pattern requiring the .wls extension.
mapfile -t TESTS < <(find tests/wolfram -maxdepth 1 -name 'test_*.wls' | sort)

if [ ${#TESTS[@]} -eq 0 ]; then
    echo "No test files found under tests/wolfram/" >&2
    exit 1
fi

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
