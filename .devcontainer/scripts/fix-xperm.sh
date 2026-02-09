#!/usr/bin/env bash
# fix-xperm.sh
# xPerm MathLink Compatibility Fix
# Patches xPerm to work without requiring GLIBC 2.38
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
XPERM_FILE="$WOLFRAM_USERBASE/Applications/xAct/xPerm/xPerm.m"

echo ""
echo "${BOLD}xPerm MathLink Compatibility Fix${RESET}"
echo "================================"
echo ""
echo "This script patches xPerm to work without requiring GLIBC 2.38."
echo "It modifies the MathLink connection logic to gracefully handle"
echo "incompatible binaries and fall back to pure Mathematica mode."
echo ""

# Check if xPerm file exists
if [[ ! -f "$XPERM_FILE" ]]; then
    print_error "xPerm.m not found at $XPERM_FILE"
    echo ""
    echo "${BOLD}Solution:${RESET}"
    echo "  Install xAct first: ${BLUE}bash .devcontainer/scripts/setup_xact.sh${RESET}"
    echo ""
    exit 1
fi

print_info "Found xPerm.m at: $XPERM_FILE"
print_info "Applying compatibility fix..."
echo ""

# Create backup if it doesn't exist
if [[ ! -f "$XPERM_FILE.backup" ]]; then
    cp "$XPERM_FILE" "$XPERM_FILE.backup"
    print_success "Created backup: xPerm.m.backup"
else
    print_info "Backup already exists: xPerm.m.backup"
fi

# Replace the problematic mathlink connection logic
# This sed command modifies the Print statement to add a FileExistsQ check
# that prevents errors when the MathLink binary is incompatible
sed -i '/Print\[If\[result,"Connection established\.","Connection failed\."\]\]/c\
    If[FileExistsQ[$xpermExecutable],\
        Check[$xpermLink=Install[$xpermExecutable],result=False];\
        Print[If[result,"Connection established.","Connection failed."]],\
        Print["MathLink binary not compatible. Using pure Mathematica mode."];\
        result=False\
    ]' "$XPERM_FILE"

echo ""
print_success "xPerm patched for compatibility"
echo ""
echo "${BOLD}What was changed:${RESET}"
echo "  • Added FileExistsQ check for MathLink executable"
echo "  • Added graceful fallback to pure Mathematica mode"
echo "  • Improved error messages for debugging"
echo ""
echo "${BOLD}Next steps:${RESET}"
echo "  • Test xPerm: ${BLUE}wolframscript -code '<<xAct\`xPerm\`; Print[\$xpermQ]'${RESET}"
echo "  • Or build MathLink: ${BLUE}bash .devcontainer/scripts/build-xperm.sh${RESET}"
echo ""