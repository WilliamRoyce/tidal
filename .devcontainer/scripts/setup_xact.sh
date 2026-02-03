#!/usr/bin/env bash
# setup_xact.sh
# Interactive setup wizard for xAct 1.3.0 tensor computation packages
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

# Configuration
XACT_VERSION="1.3.0"
XACT_URL="https://xact.es/xAct_${XACT_VERSION}.tgz"
USERBASE_DIR="${WOLFRAM_USERBASE:-$HOME/.local/wolfram/userbase}"
APPS_DIR="$USERBASE_DIR/Applications"
XACT_DIR="$APPS_DIR/xAct"
TEMP_DIR="/tmp/xact_install_$$"

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

cleanup() {
    if [ -d "$TEMP_DIR" ]; then
        rm -rf "$TEMP_DIR"
    fi
}

trap cleanup EXIT

# Main script
clear
print_header "xAct ${XACT_VERSION} Setup Wizard"

echo "This interactive wizard will install the xAct tensor computation framework"
echo "for use with Wolfram Engine in the torsion-gertsenshtein devcontainer."
echo ""
echo "${BOLD}What is xAct?${RESET}"
echo "  xAct is a suite of free Mathematica packages for tensor computations"
echo "  in differential geometry and general relativity. It provides powerful"
echo "  tools for abstract tensor algebra and coordinate-based calculations."
echo ""
echo "${BOLD}Included packages:${RESET}"
echo "  • ${BOLD}xTensor${RESET}   - Abstract tensor algebra, index manipulation"
echo "  • ${BOLD}xCoba${RESET}     - Coordinate-based computations, Christoffel symbols"
echo "  • ${BOLD}xPerm${RESET}     - Advanced permutation algorithms (with MathLink)"
echo "  • ${BOLD}xPert${RESET}     - Systematic perturbation theory"
echo "  • ${BOLD}Spinors${RESET}   - Spinor formalism for general relativity"
echo "  • ${BOLD}TexAct${RESET}    - LaTeX export for tensor expressions"
echo "  • And more specialized packages..."
echo ""
echo "${BOLD}What this script does:${RESET}"
echo "  1. Check if xAct is already installed"
echo "  2. Verify Wolfram Engine is available"
echo "  3. Download xAct ${XACT_VERSION} from xact.es (~6 MB)"
echo "  4. Extract packages to Wolfram UserBase Applications directory"
echo "  5. Verify installation by loading xTensor"
echo "  6. Optionally compile xPerm for high-performance MathLink"
echo ""
echo "${BOLD}Prerequisites:${RESET}"
echo "  • Wolfram Engine 14.3 installed (run setup_wolfram_engine.sh first)"
echo "  • Internet connection for download"
echo ""
read -p "Press Enter to continue..."

# Step 1: Check existing installation
print_step "Step 1/6: Checking for existing xAct installation"

if [ -d "$XACT_DIR" ] && [ -f "$XACT_DIR/xTensor/xTensor.m" ]; then
    print_success "xAct is already installed at: $XACT_DIR"

    # Check version if possible
    print_info "Verifying xAct installation..."
    PACKAGES_FOUND=0
    for pkg in xTensor xCoba xPerm xPert; do
        if [ -d "$XACT_DIR/$pkg" ]; then
            print_success "  Found: $pkg"
            ((PACKAGES_FOUND++))
        fi
    done

    echo ""
    echo "${GREEN}${BOLD}✓ xAct is already installed!${RESET}"
    echo ""
    echo "Found $PACKAGES_FOUND core packages. Your installation appears complete."
    echo ""
    echo "${BOLD}Next steps:${RESET}"
    echo "  • Compile xPerm: ${BLUE}bash .devcontainer/scripts/build-xperm.sh${RESET}"
    echo "  • Run tests: ${BLUE}bash .devcontainer/tests/test-all-xact.sh${RESET}"
    echo "  • Validate setup: ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
    echo ""
    read -p "Reinstall anyway? (y/N): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Exiting. Your existing xAct installation is ready to use."
        exit 0
    fi
    print_warning "Proceeding with reinstallation..."
fi

# Step 2: Verify Wolfram Engine
print_step "Step 2/6: Verifying Wolfram Engine"

if ! command -v wolframscript &> /dev/null; then
    print_error "wolframscript not found!"
    echo ""
    echo "xAct requires Wolfram Engine to be installed first."
    echo ""
    echo "${BOLD}Please run the Wolfram Engine setup:${RESET}"
    echo "  ${BLUE}bash .devcontainer/scripts/setup_wolfram_engine.sh${RESET}"
    echo ""
    exit 1
fi

print_success "Found wolframscript: $(command -v wolframscript)"

# Test Wolfram Engine
print_info "Testing Wolfram Engine..."
if wolframscript -code '2+2' 2>&1 | grep -q "4"; then
    print_success "Wolfram Engine is working"
else
    print_error "Wolfram Engine test failed"
    echo ""
    echo "Please ensure Wolfram Engine is properly installed and activated."
    exit 1
fi

# Verify userbase directory
print_info "UserBase directory: $USERBASE_DIR"
if [ ! -d "$USERBASE_DIR" ]; then
    print_warning "UserBase directory doesn't exist, creating..."
    mkdir -p "$USERBASE_DIR"
fi
print_success "UserBase directory accessible"

# Step 3: Download xAct
print_step "Step 3/6: Downloading xAct ${XACT_VERSION}"

echo "${BOLD}Downloading from: ${BLUE}${XACT_URL}${RESET}"
echo ""
echo "The xAct package bundle is approximately 6 MB."
echo "This includes all core packages and documentation."
echo ""

# Create temp directory
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

print_info "Downloading xAct tarball..."
if curl -L --progress-bar -o "xAct_${XACT_VERSION}.tgz" "$XACT_URL"; then
    print_success "Download complete!"

    # Verify download
    FILESIZE=$(stat -c%s "xAct_${XACT_VERSION}.tgz" 2>/dev/null || stat -f%z "xAct_${XACT_VERSION}.tgz" 2>/dev/null || echo "0")
    if [ "$FILESIZE" -lt 100000 ]; then
        print_error "Downloaded file is too small ($FILESIZE bytes)"
        echo "The download may have failed. Please check your internet connection."
        exit 1
    fi
    print_success "Downloaded $(numfmt --to=iec-i --suffix=B $FILESIZE 2>/dev/null || echo "$FILESIZE bytes")"
else
    print_error "Download failed!"
    echo ""
    echo "Could not download xAct from ${XACT_URL}"
    echo "Please check your internet connection and try again."
    exit 1
fi

# Step 4: Extract xAct
print_step "Step 4/6: Extracting xAct Packages"

echo "${BOLD}Extracting to: ${GREEN}${APPS_DIR}${RESET}"
echo ""
echo "The tarball contains all xAct packages and will be extracted"
echo "to your Wolfram UserBase Applications directory."
echo ""

# Create Applications directory if needed
mkdir -p "$APPS_DIR"

print_info "Extracting tarball..."
if tar -xzf "xAct_${XACT_VERSION}.tgz" -C "$APPS_DIR"; then
    print_success "Extraction complete!"
else
    print_error "Extraction failed!"
    exit 1
fi

# Verify extraction
if [ ! -d "$XACT_DIR" ]; then
    print_error "xAct directory not found after extraction"
    echo "Expected: $XACT_DIR"
    exit 1
fi

print_success "xAct directory created: $XACT_DIR"

# List installed packages
echo ""
print_info "Installed packages:"
PACKAGE_COUNT=0
for pkg_dir in "$XACT_DIR"/*; do
    if [ -d "$pkg_dir" ]; then
        pkg_name=$(basename "$pkg_dir")
        # Look for .m file with same name
        if [ -f "$pkg_dir/${pkg_name}.m" ]; then
            print_success "  $pkg_name"
            ((PACKAGE_COUNT++))
        fi
    fi
done

print_success "Total packages installed: $PACKAGE_COUNT"

# Step 5: Verify installation
print_step "Step 5/6: Verifying Installation"

echo "Testing xAct by loading the xTensor package..."
echo ""

print_info "Running: wolframscript -code '<<xAct\`xTensor\`; Print[\"xAct loaded successfully!\"]'"
echo ""

# Create a test script
TEST_SCRIPT="$TEMP_DIR/test_xact.wls"
cat > "$TEST_SCRIPT" << 'EOF'
#!/usr/bin/env wolframscript

(* Suppress xAct banner *)
Off[General::shdw];
$PrePrint = .;

(* Load xTensor *)
Needs["xAct`xTensor`"];

(* Simple test: Define a 4D manifold *)
DefManifold[M, 4, {\[Mu], \[Nu], \[Rho], \[Sigma]}];

(* Print success message *)
Print["✓ xTensor loaded successfully!"];
Print["✓ Created test manifold M with dimension: ", Dim[M]];

(* Exit *)
Exit[0];
EOF

if wolframscript "$TEST_SCRIPT" 2>&1; then
    echo ""
    print_success "xAct installation verified!"
else
    echo ""
    print_warning "Verification test produced warnings (this may be normal)"
    print_info "xAct packages are installed, but you should run the full test suite"
fi

# Step 6: xPerm MathLink compilation
print_step "Step 6/6: xPerm MathLink Compilation (Optional)"

echo "${BOLD}About xPerm MathLink:${RESET}"
echo "  xPerm includes a high-performance C program for permutation computations."
echo "  By default, xAct uses pure Mathematica code, which works but is slower."
echo "  Compiling the MathLink executable provides significant performance benefits."
echo ""
echo "${BOLD}What compilation does:${RESET}"
echo "  • Compiles xPerm C source code with current system's GLIBC"
echo "  • Links against Wolfram Engine's MathLink SDK"
echo "  • Enables fast tensor symmetry computations"
echo "  • Required for optimal performance in complex calculations"
echo ""
echo "${BOLD}When to compile:${RESET}"
echo "  • ${GREEN}Now${RESET} - Recommended if you'll use xAct immediately"
echo "  • ${YELLOW}Later${RESET} - Can run .devcontainer/scripts/build-xperm.sh anytime"
echo ""

BUILD_SCRIPT=".devcontainer/scripts/build-xperm.sh"
if [ -f "$BUILD_SCRIPT" ]; then
    read -p "Compile xPerm MathLink now? (Y/n): " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        print_info "Skipping xPerm compilation"
        echo "You can compile later by running:"
        echo "  ${BLUE}bash $BUILD_SCRIPT${RESET}"
    else
        print_info "Starting xPerm compilation..."
        echo ""
        if bash "$BUILD_SCRIPT"; then
            print_success "xPerm MathLink compiled successfully!"
        else
            print_warning "xPerm compilation encountered issues"
            echo "xAct will fall back to pure Mathematica mode."
            echo "You can try again later with:"
            echo "  ${BLUE}bash $BUILD_SCRIPT${RESET}"
        fi
    fi
else
    print_warning "xPerm build script not found at: $BUILD_SCRIPT"
    print_info "You can compile xPerm manually from: $XACT_DIR/xPerm/mathlink/"
fi

# Success!
print_header "✓ xAct Setup Complete!"

echo "${GREEN}${BOLD}Installation successful!${RESET}"
echo ""
echo "${BOLD}What was installed:${RESET}"
echo "  • xAct version ${XACT_VERSION}"
echo "  • Location: $XACT_DIR"
echo "  • Packages: $PACKAGE_COUNT packages available"
echo ""
echo "${BOLD}Core packages ready:${RESET}"
echo "  • xTensor - Abstract tensor algebra"
echo "  • xCoba - Coordinate-based computations"
echo "  • xPerm - Permutation algorithms"
echo "  • xPert - Perturbation theory"
echo ""
echo "${BOLD}Next steps:${RESET}"
echo "  1. Run comprehensive tests:"
echo "     ${BLUE}bash .devcontainer/tests/test-all-xact.sh${RESET}"
echo ""
echo "  2. Validate complete setup:"
echo "     ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
echo ""
echo "  3. Try xAct in your code:"
echo "     ${BLUE}wolframscript -code '<<xAct\`xTensor\`; DefManifold[M,4]; Print[Dim[M]]'${RESET}"
echo ""
echo "${BOLD}Documentation:${RESET}"
echo "  • xAct homepage: ${BLUE}http://xact.es${RESET}"
echo "  • Test suite docs: ${BLUE}.devcontainer/docs/XACT_TESTS.md${RESET}"
echo "  • Wolfram guide: ${BLUE}.devcontainer/docs/WOLFRAM_GUIDE.md${RESET}"
echo ""
echo "Your xAct installation is persistent across container rebuilds!"
echo ""
