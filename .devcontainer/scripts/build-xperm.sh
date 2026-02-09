#!/usr/bin/env bash
# build-xperm.sh
# xPerm MathLink Compiler and Installer
# Compiles xPerm from source with current GLIBC and installs with proper library path
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

print_warning() {
    echo "${YELLOW}⚠${RESET} $1"
}

print_error() {
    echo "${RED}✗${RESET} $1"
}

print_info() {
    echo "${BLUE}ℹ${RESET} $1"
}

print_step() {
    echo ""
    echo "${BOLD}=== $1 ===${RESET}"
    echo ""
}

# Configuration
XPERM_DIR="$WOLFRAM_USERBASE/Applications/xAct/xPerm/mathlink"
MATHLINK_DIR="/home/vscode/.local/wolfram/engine/14.3/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions"

clear
print_header "xPerm MathLink Compiler and Installer"

echo "This script compiles the xPerm MathLink binary from source code,"
echo "enabling high-performance permutation computations for tensor algebra."
echo ""
echo "${BOLD}What this script does:${RESET}"
echo "  1. Verify build prerequisites (gcc, MathLink SDK)"
echo "  2. Create backups of original files"
echo "  3. Generate MathLink template code with mprep"
echo "  4. Compile xPerm C source against current GLIBC"
echo "  5. Create wrapper script for proper library paths"
echo "  6. Test the compiled MathLink connection"
echo ""
echo "${BOLD}Why compile from source?${RESET}"
echo "  • The pre-compiled xPerm binary may require newer GLIBC"
echo "  • Compiling with system GLIBC ensures compatibility"
echo "  • Provides significant performance boost for complex calculations"
echo ""

# Step 1: Check prerequisites
print_step "Step 1/6: Checking Prerequisites"

if ! command -v gcc &> /dev/null; then
    print_warning "gcc not found. Installing build tools..."
    sudo apt update && sudo apt install -y build-essential uuid-dev
    print_success "Build tools installed"
else
    print_success "gcc compiler found: $(gcc --version | head -1)"
fi

if [[ ! -f "$MATHLINK_DIR/mprep" ]]; then
    print_error "Wolfram MathLink SDK not found at $MATHLINK_DIR"
    echo ""
    echo "${BOLD}Solution:${RESET}"
    echo "  Install Wolfram Engine 14.3: ${BLUE}bash .devcontainer/scripts/setup_wolfram_engine.sh${RESET}"
    echo ""
    exit 1
fi

print_success "MathLink SDK found: $MATHLINK_DIR"

if [[ ! -d "$XPERM_DIR" ]]; then
    print_error "xPerm package not found at $XPERM_DIR"
    echo ""
    echo "${BOLD}Solution:${RESET}"
    echo "  Install xAct: ${BLUE}bash .devcontainer/scripts/setup_xact.sh${RESET}"
    echo ""
    exit 1
fi

cd "$XPERM_DIR"
print_success "All prerequisites satisfied"
echo ""

# Step 2: Backup original files
print_step "Step 2/6: Creating Backups"

echo "Preserving original files before compilation..."
echo ""

if [[ ! -f "xperm.c.original" ]]; then
    cp xperm.c xperm.c.original
    print_success "Created: xperm.c.original"
else
    print_info "Backup already exists: xperm.c.original"
fi

if [[ -f "xperm.linux.64-bit" && ! -f "xperm.linux.64-bit.factory" ]]; then
    cp xperm.linux.64-bit xperm.linux.64-bit.factory
    print_success "Created: xperm.linux.64-bit.factory"
elif [[ -f "xperm.linux.64-bit.factory" ]]; then
    print_info "Backup already exists: xperm.linux.64-bit.factory"
fi

# Step 3: Generate MathLink template
print_step "Step 3/6: Generating MathLink Template"

echo "Running mprep to process xperm.tm template file..."
echo ""

if "$MATHLINK_DIR/mprep" xperm.tm -o xpermp.c; then
    print_success "MathLink template generated: xpermp.c"
else
    print_error "Template generation failed"
    exit 1
fi

# Step 4: Compile xPerm
print_step "Step 4/6: Compiling xPerm"

echo "Compiling C source with current system GLIBC..."
echo ""
print_info "Compiler: $(gcc --version | head -1)"
print_info "GLIBC: $(ldd --version | head -1)"
echo ""

if gcc -I"$MATHLINK_DIR" -L"$MATHLINK_DIR" -O2 xpermp.c \
    -lML64i4 -lpthread -lrt -lstdc++ -ldl -luuid -lm \
    -o xperm.linux.64-bit.compiled; then
    print_success "Compilation successful: xperm.linux.64-bit.compiled"
else
    print_error "Compilation failed"
    exit 1
fi

# Step 5: Create wrapper script
print_step "Step 5/6: Creating MathLink Wrapper"

echo "Creating wrapper script to set LD_LIBRARY_PATH for MathLink libraries..."
echo ""

cat > xperm.linux.64-bit.wrapper << 'EOF'
#!/bin/bash
# xPerm MathLink Wrapper - Auto-generated
# This wrapper sets the correct library path for MathLink shared libraries
MATHLINK_DIR="/home/vscode/.local/wolfram/engine/14.3/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions"
export LD_LIBRARY_PATH="$MATHLINK_DIR:$LD_LIBRARY_PATH"
exec "$(dirname "$0")/xperm.linux.64-bit.compiled" "$@"
EOF

chmod +x xperm.linux.64-bit.wrapper
print_success "Wrapper script created and made executable"

# Install new version
print_info "Installing new xPerm binary..."
if [[ -f "xperm.linux.64-bit" ]]; then
    mv xperm.linux.64-bit xperm.linux.64-bit.old
    print_info "Moved old binary to: xperm.linux.64-bit.old"
fi
mv xperm.linux.64-bit.wrapper xperm.linux.64-bit
print_success "Installed: xperm.linux.64-bit"

# Step 6: Test installation
print_step "Step 6/6: Testing MathLink Connection"

echo "Running integration test with Wolfram Engine..."
echo ""

if timeout 10 wolframscript -code '<<xAct`xPerm`; Print["Connection test: ", $xpermQ];' 2>/dev/null | grep -q "Connection established"; then
    echo ""
    print_header "✓ xPerm MathLink Build Complete!"

    echo "${GREEN}${BOLD}SUCCESS: xPerm advanced algorithms now available!${RESET}"
    echo ""
    echo "${BOLD}Performance features enabled:${RESET}"
    echo "  • Fast permutation group computations"
    echo "  • Strong generating sets"
    echo "  • Stabilizer chain algorithms"
    echo "  • All xPerm functions at full performance"
    echo ""
else
    echo ""
    print_warning "Installation completed but test failed"
    echo ""
    echo "${BOLD}This may be normal if:${RESET}"
    echo "  • Wolfram Engine is not activated"
    echo "  • xAct packages need to be reloaded"
    echo ""
    echo "${BOLD}Manual verification:${RESET}"
    echo "  ${BLUE}wolframscript -code '<<xAct\`xPerm\`; Print[\$xpermQ]'${RESET}"
    echo ""
fi

echo "${BOLD}Files created/modified:${RESET}"
echo "  • xperm.linux.64-bit (wrapper script)"
echo "  • xperm.linux.64-bit.compiled (compiled binary)"
echo "  • xperm.c.original (backup)"
echo "  • xperm.linux.64-bit.factory (factory backup)"
echo ""
echo "${BOLD}To restore factory version:${RESET}"
echo "  ${BLUE}cd $XPERM_DIR${RESET}"
echo "  ${BLUE}mv xperm.linux.64-bit.factory xperm.linux.64-bit${RESET}"
echo ""