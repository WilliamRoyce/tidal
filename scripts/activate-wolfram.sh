#!/usr/bin/env bash
# ==============================================================================
# Wolfram Engine Activation Script
# ==============================================================================
# This script helps activate Wolfram Engine with various methods.
#
# Usage:
#   ./scripts/activate-wolfram.sh [--interactive | --env | --check]
#
# Options:
#   --interactive  Run interactive activation (default)
#   --env          Activate using WOLFRAM_ID and WOLFRAM_PASSWORD environment variables
#   --check        Check current activation status
#
# Environment Variables (for --env mode):
#   WOLFRAM_ID        Your Wolfram ID (email)
#   WOLFRAM_PASSWORD  Your Wolfram password (use with caution)
#
# ==============================================================================

set -euo pipefail

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

# Check if Wolfram Engine is installed
check_installation() {
    if ! command -v wolframscript &> /dev/null; then
        log_error "Wolfram Engine is not installed or not in PATH"
        log_error "Run: sudo ./scripts/install-wolfram-engine.sh"
        exit 1
    fi
}

# Check current activation status
check_activation() {
    log_info "Checking activation status..."
    
    # Try to run a simple computation
    local result
    if result=$(wolframscript -code '1+1' 2>&1); then
        if [[ "$result" == "2" ]]; then
            log_info "✓ Wolfram Engine is activated and working!"
            
            # Get version info
            local version
            version=$(wolframscript -code '$VersionNumber' 2>/dev/null || echo "unknown")
            log_info "  Version: ${version}"
            
            # Get license info
            local license_type
            license_type=$(wolframscript -code '$LicenseType' 2>/dev/null || echo "unknown")
            log_info "  License Type: ${license_type}"
            
            return 0
        fi
    fi
    
    log_warn "✗ Wolfram Engine is not activated"
    log_warn "  Error: ${result}"
    return 1
}

# Interactive activation
activate_interactive() {
    log_info "Starting interactive activation..."
    log_info ""
    log_step "You will be prompted to enter your Wolfram ID credentials."
    log_step "If you don't have a Wolfram ID, create one at: https://account.wolfram.com/"
    log_info ""
    
    wolframscript -activate
    
    # Verify activation
    if check_activation; then
        log_info "Activation successful!"
        return 0
    else
        log_error "Activation may have failed. Please try again."
        return 1
    fi
}

# Activate using environment variables
activate_with_env() {
    if [[ -z "${WOLFRAM_ID:-}" ]]; then
        log_error "WOLFRAM_ID environment variable not set"
        exit 1
    fi
    
    if [[ -z "${WOLFRAM_PASSWORD:-}" ]]; then
        log_error "WOLFRAM_PASSWORD environment variable not set"
        exit 1
    fi
    
    log_info "Activating with environment credentials..."
    log_warn "Note: This method stores credentials in the command history"
    
    # Use expect-style activation if available, otherwise fall back to interactive
    # Unfortunately, wolframscript -activate doesn't accept credentials via stdin easily
    # This is a limitation of the Wolfram activation process
    
    log_warn "Environment-based activation is limited."
    log_warn "For automated activation, consider using a license server or"
    log_warn "pre-activating and copying the licensing directory."
    
    wolframscript -activate
}

# Show help for activation
show_activation_help() {
    log_info "=========================================="
    log_info "Wolfram Engine Activation Guide"
    log_info "=========================================="
    echo ""
    echo "OPTION 1: Interactive Activation (Recommended)"
    echo "  Run: wolframscript -activate"
    echo "  - Enter your Wolfram ID (email)"
    echo "  - Enter your password"
    echo ""
    echo "OPTION 2: Web-based Activation"
    echo "  1. Run: wolframscript -activate"
    echo "  2. Choose web-based activation option"
    echo "  3. Visit the URL provided"
    echo "  4. Complete activation in browser"
    echo ""
    echo "OPTION 3: License File (Enterprise)"
    echo "  - Copy mathpass file to ~/.WolframEngine/Licensing/"
    echo "  - Contact your Wolfram administrator for details"
    echo ""
    echo "CREATING A WOLFRAM ID:"
    echo "  1. Visit: https://account.wolfram.com/login/create"
    echo "  2. Create a free account"
    echo "  3. Wolfram Engine free license includes:"
    echo "     - 2GB memory limit"
    echo "     - Single machine use"
    echo "     - Non-commercial/educational use"
    echo ""
    log_info "=========================================="
}

# Main function
main() {
    local mode="interactive"
    
    # Parse arguments
    for arg in "$@"; do
        case $arg in
            --interactive|-i)
                mode="interactive"
                shift
                ;;
            --env|-e)
                mode="env"
                shift
                ;;
            --check|-c)
                mode="check"
                shift
                ;;
            --help|-h)
                show_activation_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $arg"
                echo "Usage: $0 [--interactive | --env | --check | --help]"
                exit 1
                ;;
        esac
    done
    
    check_installation
    
    case $mode in
        check)
            check_activation
            ;;
        interactive)
            activate_interactive
            ;;
        env)
            activate_with_env
            ;;
    esac
}

# Run main function
main "$@"
