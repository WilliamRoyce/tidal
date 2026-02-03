#!/bin/bash
# Master Test Runner for xAct Package Suite
# Runs all individual package tests and integration test

set -e

SCRIPT_DIR="$(dirname "$0")"
echo "=== xAct Package Suite Test Runner ==="
echo "Location: $SCRIPT_DIR"
echo ""

# Make all test scripts executable
echo "Making test scripts executable..."
chmod +x "$SCRIPT_DIR"/test-*.wls
echo "✅ Test scripts ready"
echo ""

# Function to run test with timeout and error handling
run_test() {
    local test_name="$1"
    local script="$2"
    local timeout_duration="$3"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "🧪 Running $test_name..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if timeout "$timeout_duration" "$script" 2>&1; then
        echo ""
        echo "✅ $test_name: PASSED"
    else
        echo ""
        echo "❌ $test_name: FAILED or TIMEOUT"
        return 1
    fi
    echo ""
}

# Test results tracking
declare -A test_results
total_tests=0
passed_tests=0

# Individual package tests
tests=(
    "xTensor Core Algebra|$SCRIPT_DIR/test-xtensor.wls|30s"
    "xCoba Coordinates|$SCRIPT_DIR/test-xcoba.wls|45s" 
    "xPerm Permutations|$SCRIPT_DIR/test-xperm.wls|30s"
    "xPert Perturbations|$SCRIPT_DIR/test-xpert.wls|30s"
    "Full Integration|$SCRIPT_DIR/test-integration.wls|60s"
)

# Run all tests
for test_spec in "${tests[@]}"; do
    IFS="|" read -r test_name script timeout <<< "$test_spec"
    total_tests=$((total_tests + 1))
    
    if run_test "$test_name" "$script" "$timeout"; then
        test_results["$test_name"]="PASSED"
        passed_tests=$((passed_tests + 1))
    else
        test_results["$test_name"]="FAILED"
    fi
done

# Summary report
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📊 TEST SUMMARY REPORT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Individual results
for test_name in "${!test_results[@]}"; do
    result="${test_results[$test_name]}"
    if [[ "$result" == "PASSED" ]]; then
        echo "✅ $test_name"
    else
        echo "❌ $test_name"
    fi
done

echo ""
echo "Statistics:"
echo "  Total Tests: $total_tests"
echo "  Passed: $passed_tests"
echo "  Failed: $((total_tests - passed_tests))"
echo "  Success Rate: $(( passed_tests * 100 / total_tests ))%"
echo ""

# Overall status
if [[ $passed_tests -eq $total_tests ]]; then
    echo "🎉 ALL TESTS PASSED - xAct suite fully functional!"
    echo ""
    echo "Your system is ready for:"
    echo "  • Advanced tensor calculations"
    echo "  • Coordinate-based computations"
    echo "  • High-performance permutation algorithms"
    echo "  • Systematic perturbation theory"
    echo "  • Complete General Relativity research"
    exit 0
else
    echo "⚠️  Some tests failed. Check individual results above."
    echo "Consider running: .devcontainer/build-xperm.sh"
    exit 1
fi