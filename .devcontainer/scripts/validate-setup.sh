#!/usr/bin/env bash
# validate-setup.sh
# Comprehensive validation script for Wolfram Engine and xAct setup
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

# Configuration
ENGINE_VERSION="14.3"
ENGINE_DIR="$HOME/.local/wolfram/engine/14.3"
USERBASE_DIR="$HOME/.local/wolfram/userbase"
KERNEL_PATH="$ENGINE_DIR/Executables/WolframKernel"
XACT_DIR="$USERBASE_DIR/Applications/xAct"

# Counters
TOTAL_CHECKS=0
PASSED_CHECKS=0
FAILED_CHECKS=0

# Verbose mode
VERBOSE=false
if [ "$1" == "-v" ] || [ "$1" == "--verbose" ]; then
    VERBOSE=true
fi

# Helper functions
print_check() {
    TOTAL_CHECKS=$((TOTAL_CHECKS + 1))
    printf "[%d/9] %-45s " "$TOTAL_CHECKS" "$1"
}

pass_check() {
    PASSED_CHECKS=$((PASSED_CHECKS + 1))
    echo "${GREEN}✅ PASS${RESET}"
    if [ "$VERBOSE" = true ] && [ -n "$1" ]; then
        echo "    ${BLUE}ℹ${RESET} $1"
    fi
}

fail_check() {
    FAILED_CHECKS=$((FAILED_CHECKS + 1))
    echo "${RED}❌ FAIL${RESET}"
    if [ -n "$1" ]; then
        echo "    ${RED}✗${RESET} $1"
    fi
    if [ -n "$2" ]; then
        echo "    ${YELLOW}⚙${RESET} Fix: $2"
    fi
}

print_header() {
    echo ""
    echo "${BOLD}${BLUE}================================================================================${RESET}"
    echo "${BOLD}${BLUE}  $1${RESET}"
    echo "${BOLD}${BLUE}================================================================================${RESET}"
    echo ""
}

print_section() {
    echo ""
    echo "${BOLD}${BLUE}$1${RESET}"
    echo ""
}

# Check functions
check_wolfram_installation() {
    print_check "Wolfram Engine 14.3 Installation"

    if [ ! -d "$ENGINE_DIR" ]; then
        fail_check "Wolfram Engine directory not found: $ENGINE_DIR" \
                   "Run: bash .devcontainer/scripts/setup_wolfram_engine.sh"
        return 1
    fi

    if [ ! -x "$KERNEL_PATH" ]; then
        fail_check "WolframKernel not executable: $KERNEL_PATH" \
                   "Check installation or run setup_wolfram_engine.sh"
        return 1
    fi

    pass_check "Found at: $ENGINE_DIR"
    return 0
}

check_wolfram_kernel() {
    print_check "WolframKernel Executable"

    if ! "$KERNEL_PATH" -noprompt -run 'Print["OK"]; Exit[]' > /dev/null 2>&1; then
        fail_check "WolframKernel failed to execute" \
                   "Check logs with: $KERNEL_PATH -noprompt -run 'Print[\$Version]; Exit[]'"
        return 1
    fi

    # Get version
    if VERSION_OUTPUT=$("$KERNEL_PATH" -noprompt -run 'Print[$VersionNumber]; Exit[]' 2>&1); then
        VERSION=$(echo "$VERSION_OUTPUT" | grep -oP '^\d+\.\d+' | head -1)
        if [ "$VERSION" == "$ENGINE_VERSION" ]; then
            pass_check "Version $VERSION verified"
        else
            fail_check "Version mismatch: found $VERSION, expected $ENGINE_VERSION" \
                       "Reinstall correct version with setup_wolfram_engine.sh"
            return 1
        fi
    else
        pass_check "Executable but version check failed"
    fi

    return 0
}

check_wolframscript() {
    print_check "WolframScript Configuration"

    if ! command -v wolframscript &> /dev/null; then
        fail_check "wolframscript command not found" \
                   "Check PATH or run: ln -s $ENGINE_DIR/Executables/wolframscript /usr/local/bin/"
        return 1
    fi

    if ! wolframscript -code '2+2' 2>&1 | grep -q "4"; then
        fail_check "wolframscript test failed (2+2 didn't return 4)" \
                   "Check configuration or activation status"
        return 1
    fi

    pass_check "Configured and functional"
    return 0
}

check_activation() {
    print_check "Activation Status"

    LICENSE_FILE="$USERBASE_DIR/Licensing/mathpass"
    if [ ! -f "$LICENSE_FILE" ]; then
        fail_check "License file not found: $LICENSE_FILE" \
                   "Activate with: wolframscript (enter Wolfram ID when prompted)"
        return 1
    fi

    # Test if activation works
    if wolframscript -code '$Version' > /dev/null 2>&1; then
        pass_check "Valid license file found"
    else
        fail_check "License file exists but activation failed" \
                   "Re-activate with: wolframscript"
        return 1
    fi

    # Check for backup
    if [ -d "$USERBASE_DIR/.activation_backup" ] || [ -d "$HOME/.cache/Wolfram" ]; then
        if [ "$VERBOSE" = true ]; then
            echo "    ${GREEN}✓${RESET} Activation backup found"
        fi
    fi

    return 0
}

check_xact_installation() {
    print_check "xAct 1.3.0 Packages"

    if [ ! -d "$XACT_DIR" ]; then
        fail_check "xAct directory not found: $XACT_DIR" \
                   "Run: bash .devcontainer/scripts/setup_xact.sh"
        return 1
    fi

    # Check for key packages
    MISSING_PACKAGES=()
    for pkg in xTensor xCoba xPerm xPert; do
        if [ ! -d "$XACT_DIR/$pkg" ]; then
            MISSING_PACKAGES+=("$pkg")
        fi
    done

    if [ ${#MISSING_PACKAGES[@]} -gt 0 ]; then
        fail_check "Missing packages: ${MISSING_PACKAGES[*]}" \
                   "Reinstall with: bash .devcontainer/scripts/setup_xact.sh"
        return 1
    fi

    # Count total packages
    PACKAGE_COUNT=$(find "$XACT_DIR" -maxdepth 2 -name "*.m" -type f | wc -l)
    pass_check "All core packages installed ($PACKAGE_COUNT .m files)"
    return 0
}

check_xperm_mathlink() {
    print_check "xPerm MathLink Compilation"

    XPERM_DIR="$XACT_DIR/xPerm"
    if [ ! -d "$XPERM_DIR" ]; then
        fail_check "xPerm package not found" \
                   "Install xAct first: bash .devcontainer/scripts/setup_xact.sh"
        return 1
    fi

    # Check for MathLink executable
    XPERM_BIN=$(find "$XPERM_DIR/mathlink" -name "xperm" -type f 2>/dev/null | head -1)
    if [ -z "$XPERM_BIN" ]; then
        fail_check "xPerm MathLink executable not found" \
                   "Compile with: bash .devcontainer/scripts/build-xperm.sh"
        return 1
    fi

    if [ ! -x "$XPERM_BIN" ]; then
        fail_check "xPerm executable not executable: $XPERM_BIN" \
                   "Run: chmod +x $XPERM_BIN or rebuild with build-xperm.sh"
        return 1
    fi

    pass_check "Compiled and executable"
    return 0
}

check_xact_tests() {
    print_check "xAct Test Suite"

    TEST_SCRIPT=".devcontainer/tests/test-all-xact.sh"
    if [ ! -f "$TEST_SCRIPT" ]; then
        fail_check "Test script not found: $TEST_SCRIPT" \
                   "Ensure you're in the project root directory"
        return 1
    fi

    # Run tests (quietly)
    TEMP_LOG=$(mktemp)
    if bash "$TEST_SCRIPT" > "$TEMP_LOG" 2>&1; then
        # Parse results
        PASSED=$(grep -oP '\d+(?= passed)' "$TEMP_LOG" | tail -1)
        TOTAL=$(grep -oP '\d+(?= total)' "$TEMP_LOG" | tail -1)

        rm -f "$TEMP_LOG"

        if [ "$PASSED" == "$TOTAL" ] && [ "$TOTAL" -gt 0 ]; then
            pass_check "All $TOTAL tests passed"
            return 0
        else
            fail_check "Only $PASSED/$TOTAL tests passed" \
                       "Debug with: bash $TEST_SCRIPT"
            return 1
        fi
    else
        rm -f "$TEMP_LOG"
        fail_check "Test suite failed to run" \
                   "Debug with: bash $TEST_SCRIPT"
        return 1
    fi
}

check_vscode_extensions() {
    print_check "VS Code Extensions"

    # Check if running in VS Code
    if [ -z "$VSCODE_IPC_HOOK_CLI" ] && [ -z "$VSCODE_GIT_IPC_HANDLE" ]; then
        pass_check "Skipped (not running in VS Code)"
        return 0
    fi

    # Check if code command is available
    if ! command -v code &> /dev/null; then
        fail_check "VS Code CLI not available" \
                   "Extensions may not be installed yet"
        return 1
    fi

    # Check for key extensions
    REQUIRED_EXTENSIONS=("ms-python.python" "charliermarsh.ruff" "Anthropic.claude-code")
    MISSING_EXTENSIONS=()

    for ext in "${REQUIRED_EXTENSIONS[@]}"; do
        if ! code --list-extensions 2>/dev/null | grep -q "$ext"; then
            MISSING_EXTENSIONS+=("$ext")
        fi
    done

    if [ ${#MISSING_EXTENSIONS[@]} -gt 0 ]; then
        fail_check "Missing extensions: ${MISSING_EXTENSIONS[*]}" \
                   "Run: bash .devcontainer/scripts/install-extensions-final.sh"
        return 1
    fi

    pass_check "Key extensions installed"
    return 0
}

check_python_wolframclient() {
    print_check "Python wolframclient Package"

    # Check if uv is available
    if ! command -v uv &> /dev/null; then
        fail_check "uv not found (Python package manager)" \
                   "Check devcontainer configuration"
        return 1
    fi

    # Check if wolframclient is installed in uv environment
    if uv pip list 2>/dev/null | grep -q "wolframclient"; then
        VERSION=$(uv pip list 2>/dev/null | grep "wolframclient" | awk '{print $2}')
        pass_check "Installed (version $VERSION)"
        return 0
    else
        fail_check "wolframclient not installed in uv environment" \
                   "Install with: uv pip install wolframclient"
        return 1
    fi
}

# Main validation
clear
print_header "Wolfram Engine Devcontainer Setup Validation"

if [ "$VERBOSE" = true ]; then
    echo "${BOLD}Running in verbose mode${RESET}"
    echo ""
fi

# Run all checks
check_wolfram_installation
check_wolfram_kernel
check_wolframscript
check_activation
check_xact_installation
check_xperm_mathlink
check_xact_tests
check_vscode_extensions
check_python_wolframclient

# Summary
print_section "Summary"

SUCCESS_RATE=0
if [ $TOTAL_CHECKS -gt 0 ]; then
    SUCCESS_RATE=$((PASSED_CHECKS * 100 / TOTAL_CHECKS))
fi

echo "${BOLD}Total Checks:${RESET} $TOTAL_CHECKS"
echo "${GREEN}${BOLD}Passed:${RESET}       $PASSED_CHECKS"
if [ $FAILED_CHECKS -gt 0 ]; then
    echo "${RED}${BOLD}Failed:${RESET}       $FAILED_CHECKS"
fi
echo "${BOLD}Success Rate:${RESET} $SUCCESS_RATE%"
echo ""

# Overall status
if [ $FAILED_CHECKS -eq 0 ]; then
    print_header "✓ Overall Setup Status: HEALTHY"
    echo "${GREEN}${BOLD}Your environment is fully configured and ready for development!${RESET}"
    echo ""
    echo "${BOLD}Next steps:${RESET}"
    echo "  • Start coding with Wolfram Engine and xAct"
    echo "  • Check documentation: ${BLUE}.devcontainer/docs/WOLFRAM_GUIDE.md${RESET}"
    echo "  • Try the examples: ${BLUE}.devcontainer/tests/test-*.wls${RESET}"
    echo ""
    exit 0
else
    print_header "⚠ Overall Setup Status: INCOMPLETE"
    echo "${YELLOW}${BOLD}Your environment has $FAILED_CHECKS issue(s) that need attention.${RESET}"
    echo ""
    echo "${BOLD}Recommended actions:${RESET}"
    echo ""

    # Provide specific remediation based on what failed
    if ! check_wolfram_installation 2>/dev/null; then
        echo "  1. ${BOLD}Install Wolfram Engine:${RESET}"
        echo "     ${BLUE}bash .devcontainer/scripts/setup_wolfram_engine.sh${RESET}"
        echo ""
    fi

    if ! check_xact_installation 2>/dev/null; then
        echo "  2. ${BOLD}Install xAct packages:${RESET}"
        echo "     ${BLUE}bash .devcontainer/scripts/setup_xact.sh${RESET}"
        echo ""
    fi

    if ! check_xperm_mathlink 2>/dev/null; then
        echo "  3. ${BOLD}Compile xPerm MathLink:${RESET}"
        echo "     ${BLUE}bash .devcontainer/scripts/build-xperm.sh${RESET}"
        echo ""
    fi

    if ! check_activation 2>/dev/null; then
        echo "  4. ${BOLD}Activate Wolfram Engine:${RESET}"
        echo "     ${BLUE}wolframscript${RESET}"
        echo "     (Enter your Wolfram ID credentials when prompted)"
        echo ""
    fi

    if ! check_python_wolframclient 2>/dev/null; then
        echo "  5. ${BOLD}Install Python client:${RESET}"
        echo "     ${BLUE}uv pip install wolframclient${RESET}"
        echo ""
    fi

    echo "${BOLD}For detailed troubleshooting:${RESET}"
    echo "  • Wolfram guide: ${BLUE}.devcontainer/docs/WOLFRAM_GUIDE.md${RESET}"
    echo "  • Test docs: ${BLUE}.devcontainer/docs/XACT_TESTS.md${RESET}"
    echo "  • Run health check: ${BLUE}bash .devcontainer/scripts/check-wolfram.sh${RESET}"
    echo ""
    echo "Run this script again after fixing issues: ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
    echo ""
    exit 1
fi
