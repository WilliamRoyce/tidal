#!/usr/bin/env bash
# ==============================================================================
# Wolfram Engine Installation Script
# ==============================================================================
# This script downloads and installs Wolfram Engine for headless Linux use.
# It is designed to be run during dev container creation.
#
# Usage:
#   ./scripts/install-wolfram-engine.sh [--skip-download]
#
# Options:
#   --skip-download  Skip downloading if installer already exists in third_party/
#
# Environment Variables:
#   WOLFRAM_VERSION  Wolfram Engine version (default: 14.3.0)
#   WOLFRAM_INSTALL_DIR  Installation directory (default: /usr/local/Wolfram/WolframEngine)
#
# ==============================================================================

set -euo pipefail

# Configuration
WOLFRAM_VERSION="${WOLFRAM_VERSION:-14.3.0}"
WOLFRAM_INSTALL_DIR="${WOLFRAM_INSTALL_DIR:-/usr/local/Wolfram/WolframEngine}"
WOLFRAM_SCRIPTS_DIR="/usr/local/bin"

# Derived paths
INSTALLER_NAME="WolframEngine_${WOLFRAM_VERSION}_LIN.sh"
DOWNLOAD_URL="https://account.wolfram.com/download/public/wolfram-engine/desktop/LINUX"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
THIRD_PARTY_DIR="${PROJECT_ROOT}/third_party"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
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

# Check if running as root or with sudo
check_privileges() {
    if [[ $EUID -ne 0 ]]; then
        log_error "This script must be run with sudo or as root"
        exit 1
    fi
}

# Check if Wolfram Engine is already installed
check_existing_installation() {
    if command -v wolframscript &> /dev/null; then
        local version
        version=$(wolframscript -code '$VersionNumber' 2>/dev/null || echo "unknown")
        log_info "Wolfram Engine already installed (version: ${version})"
        return 0
    fi
    return 1
}

# Install system dependencies
install_dependencies() {
    log_info "Installing system dependencies..."
    apt-get update -qq
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        xz-utils \
        unzip \
        libglib2.0-0 \
        libx11-6 \
        libxext6 \
        libxrender1 \
        libsm6 \
        libfontconfig1 \
        libasound2 \
        libfreetype6 \
        libgl1 \
        libglu1-mesa
    log_info "Dependencies installed successfully"
}

# Download Wolfram Engine installer
download_installer() {
    local installer_path="${THIRD_PARTY_DIR}/${INSTALLER_NAME}"
    
    # Create third_party directory if it doesn't exist
    mkdir -p "$THIRD_PARTY_DIR"
    
    if [[ -f "$installer_path" ]]; then
        log_info "Installer already exists at ${installer_path}"
        return 0
    fi
    
    log_info "Downloading Wolfram Engine ${WOLFRAM_VERSION}..."
    log_warn "Note: You may need to download manually from https://www.wolfram.com/engine/"
    log_warn "The direct download URL requires authentication."
    
    # Try to download (this may fail if authentication is required)
    if curl -L -o "$installer_path" "$DOWNLOAD_URL" 2>/dev/null; then
        # Check if we got a valid shell script (not an HTML error page)
        if head -1 "$installer_path" | grep -q "^#!"; then
            log_info "Download completed successfully"
            chmod +x "$installer_path"
            return 0
        else
            log_warn "Downloaded file does not appear to be the installer"
            rm -f "$installer_path"
        fi
    fi
    
    log_error "Automatic download failed. Please download manually:"
    log_error "  1. Visit https://www.wolfram.com/engine/"
    log_error "  2. Download the Linux installer"
    log_error "  3. Place it at: ${installer_path}"
    return 1
}

# Run the Wolfram Engine installer
run_installer() {
    local installer_path="${THIRD_PARTY_DIR}/${INSTALLER_NAME}"
    
    if [[ ! -f "$installer_path" ]]; then
        log_error "Installer not found at ${installer_path}"
        return 1
    fi
    
    log_info "Running Wolfram Engine installer..."
    
    # Run installer in automatic mode
    # The installer accepts these during automatic installation:
    # - Installation directory
    # - Scripts directory
    # - Whether to create symbolic links
    
    # Create installation directory
    mkdir -p "$WOLFRAM_INSTALL_DIR"
    
    # Run with auto mode and predefined answers
    # Note: The installer prompts are:
    # 1. Installation directory (default: /usr/local/Wolfram/WolframEngine/14.3)
    # 2. Scripts directory (default: /usr/local/bin)
    # We pipe these answers to stdin
    echo -e "\n\n" | bash "$installer_path" -- \
        -auto \
        -targetdir="${WOLFRAM_INSTALL_DIR}/${WOLFRAM_VERSION}" \
        -execdir="${WOLFRAM_SCRIPTS_DIR}" \
        2>/dev/null || {
            # If -auto flag doesn't work, try with heredoc for interactive prompts
            log_info "Trying interactive installation with preset answers..."
            bash "$installer_path" <<EOF


EOF
        }
    
    log_info "Wolfram Engine installation completed"
}

# Verify installation
verify_installation() {
    log_info "Verifying installation..."
    
    if ! command -v wolframscript &> /dev/null; then
        log_error "wolframscript not found in PATH"
        return 1
    fi
    
    log_info "Wolfram Engine installed successfully!"
    log_info "Location: $(command -v wolframscript)"
    
    # Check if activated (this will fail if not activated, which is expected)
    log_info ""
    log_info "=========================================="
    log_info "IMPORTANT: License Activation Required"
    log_info "=========================================="
    log_info "Run the following command to activate:"
    log_info "  wolframscript -activate"
    log_info ""
    log_info "You will need a Wolfram ID (free at wolfram.com)"
    log_info "=========================================="
    
    return 0
}

# Main installation function
main() {
    local skip_download=false
    
    # Parse arguments
    for arg in "$@"; do
        case $arg in
            --skip-download)
                skip_download=true
                shift
                ;;
            --help|-h)
                echo "Usage: $0 [--skip-download]"
                echo ""
                echo "Options:"
                echo "  --skip-download  Skip downloading if installer already exists"
                echo ""
                exit 0
                ;;
            *)
                log_error "Unknown option: $arg"
                exit 1
                ;;
        esac
    done
    
    log_info "=========================================="
    log_info "Wolfram Engine Installation Script"
    log_info "Version: ${WOLFRAM_VERSION}"
    log_info "=========================================="
    
    # Check if already installed
    if check_existing_installation; then
        log_info "Skipping installation (already installed)"
        exit 0
    fi
    
    check_privileges
    install_dependencies
    
    if [[ "$skip_download" == false ]]; then
        download_installer || {
            log_warn "Download failed, checking if installer exists..."
        }
    fi
    
    # Check if installer exists before proceeding
    if [[ ! -f "${THIRD_PARTY_DIR}/${INSTALLER_NAME}" ]]; then
        log_error "Installer not found. Please download manually and place in third_party/"
        exit 1
    fi
    
    run_installer
    verify_installation
    
    log_info "Installation complete!"
}

# Run main function
main "$@"
