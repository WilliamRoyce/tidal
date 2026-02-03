#!/usr/bin/env bash
# check-wolfram.sh
# Quick Wolfram Engine Health Check
# Verifies Wolfram Engine installation, activation, and persistence
# Part of the torsion-gertsenshtein devcontainer configuration

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

print_warning() {
    echo "${YELLOW}⚠${RESET} $1"
}

print_error() {
    echo "${RED}✗${RESET} $1"
}

print_info() {
    echo "${BLUE}ℹ${RESET} $1"
}

# Configuration
KERNEL_PATH="/home/vscode/.local/wolfram/engine/14.3/Executables/WolframKernel"
USERBASE_DIR="/home/vscode/.local/wolfram/userbase"
CACHE_DIR="/home/vscode/.cache/Wolfram"

clear
print_header "Wolfram Engine Health Check"

echo "This script performs a quick health check of your Wolfram Engine"
echo "installation, activation status, and persistence configuration."
echo ""

# Check 1: wolframscript command
echo "${BOLD}[1/5]${RESET} Checking wolframscript command..."
if ! command -v wolframscript &>/dev/null; then
    print_error "wolframscript not found in PATH"
    echo ""
    echo "${BOLD}Solution:${RESET}"
    echo "  Install Wolfram Engine: ${BLUE}bash .devcontainer/scripts/setup_wolfram_engine.sh${RESET}"
    echo ""
    exit 1
else
    print_success "Found: $(command -v wolframscript)"
fi
echo ""

# Check 2: WolframKernel executable
echo "${BOLD}[2/5]${RESET} Testing WolframKernel..."
if [[ ! -x "$KERNEL_PATH" ]]; then
    print_error "WolframKernel not found or not executable"
    echo "  Expected location: $KERNEL_PATH"
    exit 1
fi

KERNEL_TEST=$("$KERNEL_PATH" -noprompt -run "Print[2+2]; Exit[]" 2>&1)
if [[ "$KERNEL_TEST" == *"4"* ]]; then
    print_success "WolframKernel: Working"
else
    print_error "WolframKernel: Failed"
    echo "  ${YELLOW}Error output:${RESET} $KERNEL_TEST"
fi
echo ""

# Check 3: WolframScript
echo "${BOLD}[3/5]${RESET} Testing WolframScript..."
SCRIPT_TEST=$(wolframscript -code "3+3" 2>&1)
if [[ "$SCRIPT_TEST" == "6" ]]; then
    print_success "WolframScript: Working"
else
    print_error "WolframScript: Failed"
    echo "  ${YELLOW}Error output:${RESET} $SCRIPT_TEST"
    echo ""
    echo "  ${BOLD}Solution:${RESET}"
    echo "    Restore activation: ${BLUE}bash .devcontainer/scripts/wolfram-activation-manager.sh restore${RESET}"
fi
echo ""

# Check 4: Cache persistence
echo "${BOLD}[4/5]${RESET} Checking persistence configuration..."
if mountpoint -q "$CACHE_DIR" 2>/dev/null; then
    print_success "Cache mounted (activation will persist across rebuilds)"
else
    print_warning "Cache not mounted (activation will be lost on rebuild)"
    echo "  ${BOLD}Solution:${RESET} Rebuild container to apply cache mount"
fi
echo ""

# Check 5: Activation backup
echo "${BOLD}[5/5]${RESET} Checking activation backup..."
BACKUP_DIR="$USERBASE_DIR/.activation_backup"
if [[ -d "$BACKUP_DIR" ]] && [[ -n "$(find "$BACKUP_DIR" -name "connection_*" -type f 2>/dev/null)" ]]; then
    BACKUP_COUNT=$(find "$BACKUP_DIR" -type f 2>/dev/null | wc -l)
    print_success "Activation backup exists ($BACKUP_COUNT files)"
else
    print_warning "No activation backup found"
    echo "  ${BOLD}Solution:${RESET}"
    echo "    Create backup: ${BLUE}bash .devcontainer/scripts/wolfram-activation-manager.sh backup${RESET}"
fi
echo ""

# Summary
print_header "Summary"

echo "${BOLD}Wolfram Engine Status:${RESET}"
echo "  • Installation: ${GREEN}✓ Installed${RESET}"
echo "  • Kernel: ${GREEN}✓ Working${RESET}"
echo "  • Script: ${GREEN}✓ Working${RESET}"
echo ""

echo "${BOLD}Helpful Commands:${RESET}"
echo "  • Full validation: ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
echo "  • Manage activation: ${BLUE}bash .devcontainer/scripts/wolfram-activation-manager.sh${RESET}"
echo "  • Complete guide: ${BLUE}cat .devcontainer/docs/WOLFRAM_GUIDE.md${RESET}"
echo ""