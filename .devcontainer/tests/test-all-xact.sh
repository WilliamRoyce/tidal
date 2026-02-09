#!/usr/bin/env bash
# test-all-xact.sh
# Master Test Runner for xAct Package Suite
# Runs all individual package tests and integration test
# Part of the TIDAL devcontainer configuration

set -e

# Colors for output
if command -v tput &> /dev/null && [ -t 1 ]; then
    GREEN=$(tput setaf 2)
    YELLOW=$(tput setaf 3)
    RED=$(tput setaf 1)
    BLUE=$(tput setaf 4)
    BOLD=$(tput bold)
    RESET=$(tput sgr0)
else
    GREEN=""
    YELLOW=""
    RED=""
    BLUE=""
    BOLD=""
    RESET=""
fi

# Helper functions
print_header() {
    echo ""
    echo "${BOLD}${BLUE}================================================================================${RESET}"
    echo "${BOLD}${BLUE}  $1${RESET}"
    echo "${BOLD}${BLUE}================================================================================${RESET}"
    echo ""
}

print_success() {
    echo "${GREEN}✓${RESET} $1"
}

print_error() {
    echo "${RED}✗${RESET} $1"
}

print_info() {
    echo "${BLUE}ℹ${RESET} $1"
}

# Configuration
SCRIPT_DIR="$(dirname "$0")"

clear
print_header "xAct Package Suite Test Runner"

echo "This script runs comprehensive tests for all xAct packages to verify"
echo "that the tensor computation framework is working correctly."
echo ""
echo "${BOLD}Tests included:${RESET}"
echo "  • xTensor - Core abstract tensor algebra"
echo "  • xCoba - Coordinate-based computations"
echo "  • xPerm - Permutation algorithms (including MathLink)"
echo "  • xPert - Perturbation theory"
echo "  • Integration test - Full workflow verification"
echo ""
echo "Location: ${BLUE}$SCRIPT_DIR${RESET}"
echo ""

# Make all test scripts executable
print_info "Preparing test scripts..."
chmod +x "$SCRIPT_DIR"/test-*.wls
print_success "Test scripts ready"
echo ""

# Function to run test with timeout and error handling
run_test() {
    local test_name="$1"
    local script="$2"
    local timeout_duration="$3"

    echo ""
    echo "${BOLD}${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
    echo "${BOLD}Running: $test_name${RESET}"
    echo "${BOLD}${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
    echo ""

    if timeout "$timeout_duration" "$script" 2>&1; then
        echo ""
        print_success "$test_name: PASSED"
    else
        echo ""
        print_error "$test_name: FAILED or TIMEOUT"
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
print_header "Test Summary Report"

# Individual results
echo "${BOLD}Individual Test Results:${RESET}"
echo ""
for test_name in "${!test_results[@]}"; do
    result="${test_results[$test_name]}"
    if [[ "$result" == "PASSED" ]]; then
        print_success "$test_name"
    else
        print_error "$test_name"
    fi
done

echo ""
echo "${BOLD}Statistics:${RESET}"
echo "  Total Tests: ${BOLD}$total_tests${RESET}"
echo "  ${GREEN}Passed: $passed_tests${RESET}"
if [[ $((total_tests - passed_tests)) -gt 0 ]]; then
    echo "  ${RED}Failed: $((total_tests - passed_tests))${RESET}"
fi
echo "  Success Rate: ${BOLD}$(( passed_tests * 100 / total_tests ))%${RESET}"
echo ""

# Overall status
if [[ $passed_tests -eq $total_tests ]]; then
    print_header "✓ All Tests Passed!"

    echo "${GREEN}${BOLD}xAct suite is fully functional!${RESET}"
    echo ""
    echo "${BOLD}Your system is ready for:${RESET}"
    echo "  • Advanced tensor calculations"
    echo "  • Coordinate-based computations"
    echo "  • High-performance permutation algorithms"
    echo "  • Systematic perturbation theory"
    echo "  • Complete General Relativity research"
    echo ""
    echo "${BOLD}Next steps:${RESET}"
    echo "  • Start using xAct in your research"
    echo "  • Check examples: ${BLUE}.devcontainer/tests/test-*.wls${RESET}"
    echo "  • Read documentation: ${BLUE}http://xact.es${RESET}"
    echo ""
    exit 0
else
    print_header "⚠ Some Tests Failed"

    echo "${YELLOW}${BOLD}$((total_tests - passed_tests)) test(s) failed. Check results above.${RESET}"
    echo ""
    echo "${BOLD}Common solutions:${RESET}"
    echo "  • Compile xPerm MathLink: ${BLUE}bash .devcontainer/scripts/build-xperm.sh${RESET}"
    echo "  • Verify Wolfram activation: ${BLUE}wolframscript -code '2+2'${RESET}"
    echo "  • Check health: ${BLUE}bash .devcontainer/scripts/check-wolfram.sh${RESET}"
    echo ""
    exit 1
fi