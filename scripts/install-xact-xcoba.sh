#!/usr/bin/env bash
# ==============================================================================
# xAct & xCoba Installation Script
# ==============================================================================
# This script installs xAct and xCoba packages for symbolic tensor algebra in
# the Wolfram Engine user applications directory.
#
# It includes GLIBC compatibility fixes by recompiling the xPerm binary.
#
# Usage:
#   ./scripts/install-xact-xcoba.sh [--version VERSION]
#
# Options:
#   --version VERSION  xAct version to install (default: 1.2.1, compatible with GLIBC 2.36)
#
# ==============================================================================

set -euo pipefail

# Configuration
XACT_VERSION="${XACT_VERSION:-1.2.1}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WOLFRAM_USER_DIR=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

# Check if Wolfram Engine is available
check_wolfram() {
    if ! command -v wolframscript &> /dev/null; then
        log_error "Wolfram Engine is not installed or not activated"
        log_error "Run: sudo ./scripts/install-wolfram-engine.sh"
        log_error "Then: ./scripts/activate-wolfram.sh"
        exit 1
    fi
    
    # Test basic functionality
    if ! wolframscript -code "1+1" >/dev/null 2>&1; then
        log_error "Wolfram Engine is not properly activated"
        log_error "Run: ./scripts/activate-wolfram.sh"
        exit 1
    fi
    
    log_info "Wolfram Engine is available and activated"
}

# Get Wolfram user base directory
get_wolfram_user_dir() {
    WOLFRAM_USER_DIR=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n' | tr -d '\r')
    
    if [[ -z "$WOLFRAM_USER_DIR" ]]; then
        log_error "Could not determine Wolfram user directory"
        exit 1
    fi
    
    log_info "Wolfram user directory: ${WOLFRAM_USER_DIR}"
}

# Install development dependencies for recompiling xPerm binary
install_build_deps() {
    log_step "Installing build dependencies for xPerm binary compilation..."
    
    if command -v apt-get &> /dev/null; then
        sudo apt-get update -qq
        sudo apt-get install -y --no-install-recommends uuid-dev gcc build-essential
    else
        log_warn "Package manager not recognized. Ensure gcc, build-essential, and uuid-dev are installed"
    fi
    
    log_info "Build dependencies installed"
}

# Download and extract xAct
download_xact() {
    local apps_dir="${WOLFRAM_USER_DIR}/Applications"
    mkdir -p "$apps_dir"
    
    cd "$apps_dir"
    
    log_step "Downloading xAct version ${XACT_VERSION}..."
    
    local archive="xAct_${XACT_VERSION}.tgz"
    local url="https://xact.es/download/${archive}"
    
    if [[ -d "xAct" ]]; then
        log_warn "Existing xAct installation found. Removing..."
        rm -rf xAct *.tgz* 2>/dev/null || true
    fi
    
    if ! wget -q --show-progress "$url"; then
        log_error "Failed to download xAct from ${url}"
        exit 1
    fi
    
    log_info "Extracting xAct..."
    if ! tar -xzf "$archive" 2>/dev/null; then
        log_error "Failed to extract ${archive}"
        exit 1
    fi
    
    log_info "xAct downloaded and extracted successfully"
}

# Recompile xPerm binary for GLIBC compatibility
recompile_xperm() {
    local mathlink_dir="${WOLFRAM_USER_DIR}/Applications/xAct/xPerm/mathlink"
    
    if [[ ! -d "$mathlink_dir" ]]; then
        log_error "xPerm mathlink directory not found: ${mathlink_dir}"
        exit 1
    fi
    
    cd "$mathlink_dir"
    
    log_step "Checking GLIBC compatibility of existing xPerm binary..."
    
    # Check if existing binary works
    if ldd xperm.linux.64-bit 2>/dev/null | grep -q "not found"; then
        log_warn "Existing binary has dependency issues. Recompiling..."
    else
        # Test if binary actually works (might still have GLIBC version issues)
        if timeout 5 ./xperm.linux.64-bit >/dev/null 2>&1; then
            log_info "Existing xPerm binary is compatible"
            return 0
        else
            log_warn "Binary compatibility test failed. Recompiling..."
        fi
    fi
    
    log_step "Finding MathLink compiler..."
    
    # Find mcc (MathLink compiler)
    local mcc_path=""
    for path in "/usr/local/Wolfram/WolframEngine/*/Executables/mcc" "/usr/local/Wolfram/Mathematica/*/Executables/mcc"; do
        if [[ -f $path ]]; then
            mcc_path="$path"
            break
        fi
    done
    
    if [[ -z "$mcc_path" ]]; then
        log_error "MathLink compiler (mcc) not found"
        exit 1
    fi
    
    log_info "Using MathLink compiler: ${mcc_path}"
    
    log_step "Recompiling xPerm binary..."
    
    # Backup original
    if [[ -f "xperm.linux.64-bit" ]]; then
        cp xperm.linux.64-bit xperm.linux.64-bit.original
    fi
    
    # Compile new binary
    if "$mcc_path" xperm.tm -luuid -O3 -o xperm.linux.64-bit.new 2>&1; then
        mv xperm.linux.64-bit.new xperm.linux.64-bit
        chmod +x xperm.linux.64-bit
        log_info "xPerm binary recompiled successfully"
    else
        log_error "Failed to recompile xPerm binary"
        exit 1
    fi
}

# Test xAct/xCoba installation
test_installation() {
    log_step "Testing xAct/xCoba installation..."
    
    local test_code='
    Quiet[Needs["xAct`xCoba`"]];
    Print["xCoba version: ", $xAct`xCobaVersionNumber];
    DefManifold[TestM, 4, IndexRange[a, z]];
    Print["✓ xAct/xCoba installation test successful"];
    '
    
    if wolframscript -code "$test_code" 2>/dev/null | grep -q "installation test successful"; then
        log_info "✅ xAct/xCoba installation verified"
        return 0
    else
        log_error "❌ xAct/xCoba installation test failed"
        return 1
    fi
}

# Show usage information
show_help() {
    log_info "=========================================="
    log_info "xAct & xCoba Installation"
    log_info "=========================================="
    echo ""
    echo "This script installs xAct and xCoba packages for symbolic"
    echo "tensor algebra computations in Wolfram Engine."
    echo ""
    echo "Features:"
    echo "  • Downloads official xAct package"
    echo "  • Recompiles xPerm binary for GLIBC compatibility"  
    echo "  • Installs to user Applications directory"
    echo "  • Verifies installation with test"
    echo ""
    echo "Usage:"
    echo "  $0 [--version VERSION] [--help]"
    echo ""
    echo "Options:"
    echo "  --version VERSION  xAct version (default: ${XACT_VERSION})"
    echo "  --help, -h         Show this help"
    echo ""
    echo "Requirements:"
    echo "  • Wolfram Engine installed and activated"
    echo "  • gcc, build-essential, uuid-dev packages"
    echo ""
    log_info "=========================================="
}

# Main installation function
main() {
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --version)
                XACT_VERSION="$2"
                shift 2
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
    
    log_info "=========================================="
    log_info "xAct & xCoba Installation Script"
    log_info "Version: ${XACT_VERSION}"
    log_info "=========================================="
    
    check_wolfram
    get_wolfram_user_dir
    install_build_deps
    download_xact
    recompile_xperm
    
    if test_installation; then
        log_info ""
        log_info "🎉 Installation completed successfully!"
        log_info ""
        log_info "Usage examples:"
        log_info '  wolframscript -code "Needs[\"xAct\`xCoba\`\"]; DefManifold[M,4]"'
        log_info "  Check examples at: ${WOLFRAM_USER_DIR}/Applications/xAct/Documentation/"
        log_info ""
    else
        log_error "Installation completed but verification failed"
        exit 1
    fi
}

# Run main function
main "$@"