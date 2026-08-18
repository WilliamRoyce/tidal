#!/usr/bin/env bash
# wolfram-activation-manager.sh
# Wolfram Activation Backup Manager
# Saves and restores Wolfram Engine activation to persistent storage
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

# Configuration
USERBASE_DIR="$HOME/.local/wolfram/userbase"
BACKUP_DIR="$USERBASE_DIR/.activation_backup"
CACHE_DIR="$HOME/.cache/Wolfram"

# Verify prerequisites
if [[ ! -d "$USERBASE_DIR" ]]; then
    clear
    print_header "Error: Wolfram Engine Not Found"
    print_error "Wolfram userbase directory not found at: $USERBASE_DIR"
    echo ""
    echo "${BOLD}This usually means:${RESET}"
    echo "  • Wolfram Engine is not installed"
    echo "  • Userbase mount is not configured correctly"
    echo ""
    echo "${BOLD}Solution:${RESET}"
    echo "  Install Wolfram Engine: ${BLUE}bash .devcontainer/scripts/setup_wolfram_engine.sh${RESET}"
    echo ""
    exit 1
fi

# Function to backup activation
backup_activation() {
    print_header "Backup Activation"

    echo "Backing up Wolfram Engine activation data to persistent storage..."
    echo ""

    if [[ ! -d "$CACHE_DIR/WolframScript" ]]; then
        print_error "No activation data found in $CACHE_DIR"
        echo ""
        echo "${BOLD}This usually means:${RESET}"
        echo "  • Wolfram Engine has not been activated yet"
        echo "  • Activation was not successful"
        echo ""
        echo "${BOLD}Solution:${RESET}"
        echo "  Activate first: ${BLUE}wolframscript${RESET}"
        echo ""
        return 1
    fi

    mkdir -p "$BACKUP_DIR"
    cp -r "$CACHE_DIR"/* "$BACKUP_DIR/"

    FILE_COUNT=$(find "$BACKUP_DIR" -type f | wc -l)
    print_success "Activation backed up to: $BACKUP_DIR"
    print_success "Backed up $FILE_COUNT files"
    echo ""

    if [ "$VERBOSE" = true ]; then
        echo "${BOLD}Backup contains:${RESET}"
        find "$BACKUP_DIR" -type f | sed 's|^|     |'
        echo ""
    fi

    print_info "This backup will survive container rebuilds!"
    return 0
}

# Function to restore activation
restore_activation() {
    print_header "Restore Activation"

    echo "Restoring Wolfram Engine activation from persistent storage..."
    echo ""

    if [[ ! -d "$BACKUP_DIR" ]] || [[ ! -d "$BACKUP_DIR/WolframScript" ]]; then
        print_error "No activation backup found at $BACKUP_DIR"
        echo ""
        echo "${BOLD}This usually means:${RESET}"
        echo "  • Backup has never been created"
        echo "  • Backup directory was deleted"
        echo ""
        echo "${BOLD}Solution:${RESET}"
        echo "  1. Activate Wolfram Engine: ${BLUE}wolframscript${RESET}"
        echo "  2. Create backup: ${BLUE}$0 backup${RESET}"
        echo ""
        return 1
    fi

    mkdir -p "$CACHE_DIR"
    cp -r "$BACKUP_DIR"/* "$CACHE_DIR/"

    FILE_COUNT=$(find "$CACHE_DIR" -type f 2>/dev/null | wc -l)
    print_success "Activation restored from backup"
    print_success "Restored $FILE_COUNT files"
    echo ""

    if [ "$VERBOSE" = true ]; then
        echo "${BOLD}Restored files:${RESET}"
        find "$CACHE_DIR" -type f | sed 's|^|     |'
        echo ""
    fi

    # Test if restoration worked
    print_info "Testing restored activation..."
    if timeout 10s wolframscript -code "2+2" 2>/dev/null | grep -q "4"; then
        print_success "Wolfram Engine is now working!"
    else
        print_warning "Restoration completed but test failed"
        echo "  You may need to re-activate: ${BLUE}wolframscript${RESET}"
    fi

    return 0
}

# Function to check status
check_status() {
    print_header "Activation Status"

    echo "Checking Wolfram Engine activation and persistence configuration..."
    echo ""

    # Check 1: Cache directory mount
    echo "${BOLD}[1/4]${RESET} Cache directory mount status:"
    if mountpoint -q "$CACHE_DIR" 2>/dev/null; then
        print_success "Mounted (will persist through rebuilds)"
    else
        print_warning "Not mounted (activation will be lost on rebuild)"
        echo "  ${BOLD}Solution:${RESET} Rebuild container to apply cache mount"
    fi
    echo ""

    # Check 2: Activation data in cache
    echo "${BOLD}[2/4]${RESET} Activation data in cache:"
    if [[ -d "$CACHE_DIR/WolframScript" ]] && [[ -n "$(find "$CACHE_DIR/WolframScript" -name "connection_*" -type f 2>/dev/null)" ]]; then
        FILE_COUNT=$(find "$CACHE_DIR/WolframScript" -name "connection_*" -type f 2>/dev/null | wc -l)
        print_success "Present ($FILE_COUNT connection files)"
    else
        print_error "Not found"
        echo "  ${BOLD}Solution:${RESET} Activate with ${BLUE}wolframscript${RESET}"
    fi
    echo ""

    # Check 3: Backup exists
    echo "${BOLD}[3/4]${RESET} Activation backup:"
    if [[ -d "$BACKUP_DIR/WolframScript" ]]; then
        FILE_COUNT=$(find "$BACKUP_DIR" -type f | wc -l)
        print_success "Exists ($FILE_COUNT files)"
        print_info "Location: $BACKUP_DIR"
    else
        print_error "Not found"
        echo "  ${BOLD}Solution:${RESET} Create backup with ${BLUE}$0 backup${RESET}"
    fi
    echo ""

    # Check 4: WolframScript works
    echo "${BOLD}[4/4]${RESET} WolframScript functionality:"
    if timeout 10s wolframscript -code "2+2" 2>/dev/null | grep -q "4"; then
        print_success "Working correctly"
    else
        print_error "Not working (needs activation)"
        echo "  ${BOLD}Solution:${RESET} Activate with ${BLUE}wolframscript${RESET}"
    fi
    echo ""
}

# Parse verbose flag
VERBOSE=false
COMMAND="${1:-status}"
if [[ "$1" == "-v" ]] || [[ "$2" == "-v" ]]; then
    VERBOSE=true
    if [[ "$1" == "-v" ]]; then
        COMMAND="${2:-status}"
    fi
fi

# Execute command
case "$COMMAND" in
    "backup")
        backup_activation
        ;;

    "restore")
        restore_activation
        ;;

    "status"|"check"|"")
        check_status
        ;;

    "help"|"-h"|"--help")
        clear
        print_header "Wolfram Activation Backup Manager - Help"

        echo "${BOLD}USAGE:${RESET}"
        echo "  $0 [COMMAND] [-v]"
        echo ""

        echo "${BOLD}COMMANDS:${RESET}"
        echo "  ${GREEN}backup${RESET}   - Save current activation to persistent storage"
        echo "  ${GREEN}restore${RESET}  - Restore activation from persistent storage"
        echo "  ${GREEN}status${RESET}   - Check current activation status (default)"
        echo "  ${GREEN}help${RESET}     - Show this help message"
        echo ""

        echo "${BOLD}OPTIONS:${RESET}"
        echo "  ${GREEN}-v${RESET}       - Verbose output (show detailed file lists)"
        echo ""

        echo "${BOLD}EXAMPLES:${RESET}"
        echo "  $0                  # Check activation status"
        echo "  $0 backup           # Back up current activation"
        echo "  $0 restore -v       # Restore with verbose output"
        echo ""

        echo "${BOLD}BACKGROUND:${RESET}"
        echo "  Wolfram activation is stored in ~/.cache/Wolfram which is normally"
        echo "  lost on container rebuild. This script backs up activation data to"
        echo "  the mounted userbase directory so it persists across rebuilds."
        echo ""

        echo "${BOLD}WORKFLOW:${RESET}"
        echo "  1. Install and activate Wolfram Engine"
        echo "  2. Run: ${BLUE}$0 backup${RESET}"
        echo "  3. On container rebuild: ${BLUE}$0 restore${RESET}"
        echo ""

        echo "${BOLD}SEE ALSO:${RESET}"
        echo "  • Health check: ${BLUE}bash .devcontainer/scripts/check-wolfram.sh${RESET}"
        echo "  • Full guide: ${BLUE}cat .devcontainer/docs/WOLFRAM_GUIDE.md${RESET}"
        echo ""
        exit 0
        ;;

    *)
        clear
        print_header "Error: Unknown Command"
        print_error "Unknown command: $COMMAND"
        echo ""
        echo "${BOLD}Valid commands:${RESET} backup, restore, status, help"
        echo ""
        echo "Run ${BLUE}$0 help${RESET} for usage information"
        echo ""
        exit 1
        ;;
esac

# Footer tips for status command
if [[ "$COMMAND" == "status" ]] || [[ "$COMMAND" == "check" ]] || [[ "$COMMAND" == "" ]]; then
    echo "${BOLD}${BLUE}───────────────────────────────────────────────────────────────────────────────${RESET}"
    echo ""
    echo "${BOLD}Quick Reference:${RESET}"
    echo "  • Backup activation: ${BLUE}$0 backup${RESET}"
    echo "  • Restore activation: ${BLUE}$0 restore${RESET}"
    echo "  • Activate Wolfram: ${BLUE}wolframscript${RESET}"
    echo "  • Health check: ${BLUE}bash .devcontainer/scripts/check-wolfram.sh${RESET}"
    echo "  • Full validation: ${BLUE}bash .devcontainer/scripts/validate-setup.sh${RESET}"
    echo ""
fi