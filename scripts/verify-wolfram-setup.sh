#!/usr/bin/env bash
# ==============================================================================
# Verification Script for Wolfram Engine and xAct/xCoba Installation
# ==============================================================================
# This script checks that all Wolfram-related components are properly installed
# and working. It should be run after container creation or to diagnose issues.
#
# Usage:
#   ./scripts/verify-wolfram-setup.sh
#
# Exit Codes:
#   0 - All checks passed
#   1 - One or more checks failed
#
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ERRORS=0
WARNINGS=0

log_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

log_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
    ((ERRORS++)) || true
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
    ((WARNINGS++)) || true
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# Check 1: Wolfram Engine binary available
check_wolfram_binary() {
    log_info "Checking Wolfram Engine installation..."
    
    if command -v wolframscript &> /dev/null; then
        local path
        path=$(command -v wolframscript)
        log_pass "wolframscript found at: ${path}"
        return 0
    else
        log_fail "wolframscript not found in PATH"
        return 1
    fi
}

# Check 2: Wolfram Engine activation
check_wolfram_activation() {
    log_info "Checking Wolfram Engine activation..."
    
    local result
    if result=$(wolframscript -code '1+1' 2>&1); then
        if [[ "$result" == "2" ]]; then
            log_pass "Wolfram Engine is activated and responding correctly"
            
            # Get version info
            local version
            version=$(wolframscript -code '$VersionNumber' 2>/dev/null || echo "unknown")
            log_info "  Version: ${version}"
            
            local license
            license=$(wolframscript -code '$LicenseType' 2>/dev/null || echo "unknown")
            log_info "  License: ${license}"
            return 0
        else
            log_fail "Wolfram Engine returned unexpected result: ${result}"
            return 1
        fi
    else
        log_fail "Wolfram Engine not activated or error occurred"
        log_info "  Run: ./scripts/activate-wolfram.sh"
        return 1
    fi
}

# Check 3: User Applications directory exists
check_applications_dir() {
    log_info "Checking Wolfram Applications directory..."
    
    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local apps_dir="${user_dir}/Applications"
    
    if [[ -d "$apps_dir" ]]; then
        log_pass "Applications directory exists: ${apps_dir}"
        return 0
    else
        log_fail "Applications directory not found: ${apps_dir}"
        return 1
    fi
}

# Check 4: xAct installation
check_xact_installed() {
    log_info "Checking xAct installation..."
    
    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local xact_dir="${user_dir}/Applications/xAct"
    
    if [[ -d "$xact_dir" ]]; then
        log_pass "xAct directory found: ${xact_dir}"
        
        # Check for key packages
        local packages=("xCore" "xPerm" "xTensor" "xCoba")
        for pkg in "${packages[@]}"; do
            if [[ -d "${xact_dir}/${pkg}" ]]; then
                log_pass "  Package ${pkg} found"
            else
                log_fail "  Package ${pkg} missing"
            fi
        done
        return 0
    else
        log_fail "xAct directory not found: ${xact_dir}"
        log_info "  Run: ./scripts/install-xact-xcoba.sh"
        return 1
    fi
}

# Check 5: xPerm binary compatibility
check_xperm_binary() {
    log_info "Checking xPerm binary compatibility..."
    
    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local xperm_binary="${user_dir}/Applications/xAct/xPerm/mathlink/xperm.linux.64-bit"
    
    if [[ ! -f "$xperm_binary" ]]; then
        log_fail "xPerm binary not found: ${xperm_binary}"
        return 1
    fi
    
    # Check for GLIBC compatibility
    local ldd_output
    if ldd_output=$(ldd "$xperm_binary" 2>&1); then
        if echo "$ldd_output" | grep -q "not found"; then
            log_fail "xPerm binary has missing dependencies"
            log_info "  Run: ./scripts/install-xact-xcoba.sh to recompile"
            return 1
        else
            log_pass "xPerm binary dependencies satisfied"
            return 0
        fi
    else
        log_warn "Could not check xPerm binary dependencies"
        return 0
    fi
}

# Check 6: Load xAct packages
check_xact_loads() {
    log_info "Checking xAct package loading..."
    
    local test_code='
    Quiet[
      Needs["xAct`xCore`"];
      Needs["xAct`xPerm`"];
      Needs["xAct`xTensor`"];
      Needs["xAct`xCoba`"];
    , {General::stop, UpSetDelayed::write, SetDelayed::write}];
    Print["PACKAGES_LOADED"];
    '
    
    local result
    if result=$(timeout 60 wolframscript -code "$test_code" 2>&1); then
        if echo "$result" | grep -q "PACKAGES_LOADED"; then
            log_pass "All xAct packages load successfully"
            
            # Check for connection to xPerm external executable
            if echo "$result" | grep -q "Connection established"; then
                log_pass "  xPerm external executable connected"
            elif echo "$result" | grep -q "GLIBC"; then
                log_fail "  xPerm has GLIBC compatibility issues"
                return 1
            fi
            return 0
        else
            log_fail "xAct packages failed to load"
            log_info "  Output: ${result}"
            return 1
        fi
    else
        log_fail "Timeout or error loading xAct packages"
        return 1
    fi
}

# Check 7: Run smoke test
check_smoke_test() {
    log_info "Running xAct/xCoba smoke test..."
    
    local smoke_script="${SCRIPT_DIR}/xact_smoke.wl"
    
    if [[ ! -f "$smoke_script" ]]; then
        log_warn "Smoke test script not found: ${smoke_script}"
        return 0
    fi
    
    local result
    if result=$(timeout 120 wolframscript -file "$smoke_script" 2>&1); then
        if echo "$result" | grep -q "SMOKE TEST PASSED"; then
            log_pass "Smoke test passed"
            return 0
        else
            log_fail "Smoke test did not complete successfully"
            log_info "  Last lines of output:"
            echo "$result" | tail -10
            return 1
        fi
    else
        log_fail "Smoke test timed out or encountered an error"
        return 1
    fi
}

# Main verification function
main() {
    echo ""
    echo "========================================"
    echo "Wolfram Engine & xAct/xCoba Verification"
    echo "========================================"
    echo ""
    
    # Run all checks
    check_wolfram_binary || true
    echo ""
    
    check_wolfram_activation || true
    echo ""
    
    check_applications_dir || true
    echo ""
    
    check_xact_installed || true
    echo ""
    
    check_xperm_binary || true
    echo ""
    
    check_xact_loads || true
    echo ""
    
    check_smoke_test || true
    echo ""
    
    # Summary
    echo "========================================"
    echo "Verification Summary"
    echo "========================================"
    
    if [[ $ERRORS -eq 0 ]]; then
        log_pass "All checks passed!"
        if [[ $WARNINGS -gt 0 ]]; then
            log_warn "${WARNINGS} warning(s) noted"
        fi
        exit 0
    else
        log_fail "${ERRORS} check(s) failed"
        if [[ $WARNINGS -gt 0 ]]; then
            log_warn "${WARNINGS} warning(s) noted"
        fi
        echo ""
        log_info "To fix issues, try:"
        log_info "  1. sudo ./scripts/install-wolfram-engine.sh"
        log_info "  2. ./scripts/activate-wolfram.sh"
        log_info "  3. ./scripts/install-xact-xcoba.sh"
        exit 1
    fi
}

# Run main function
main "$@"
