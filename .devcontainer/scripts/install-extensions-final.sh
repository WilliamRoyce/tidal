#!/usr/bin/env bash
# install-extensions-final.sh
# VS Code Extension Installer
# Installs all required VS Code extensions after server initialization
# Part of the TIDAL devcontainer configuration

set -e

# Derive the repo root -- never assume a particular clone location.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

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
LOG_FILE="$HOME/.vscode-server/data/logs/install-extensions.log"
mkdir -p "$(dirname "$LOG_FILE")"

# Also output to both terminal and log file
exec > >(tee -a "$LOG_FILE") 2>&1

clear
print_header "VS Code Extension Installer"

echo "Started: ${BOLD}$(date)${RESET}"
echo ""
echo "This script installs all required VS Code extensions for the"
echo "tidal development environment."
echo ""

# Extension list with categories
EXTENSIONS=(
    "ms-python.python"
    "ms-python.debugpy"
    "esbenp.prettier-vscode"
    "streetsidesoftware.code-spell-checker"
    "tamasfe.even-better-toml"
    "charliermarsh.ruff"
    "peakchen90.open-html-in-browser"
    "github.vscode-github-actions"
    "GitHub.copilot"
    "GitHub.copilot-chat"
    "github.vscode-pull-request-github"
    "Anthropic.claude-code"
    "WolframResearch.wolfram"
    "lsp-wl.lsp-wl-client"
)

echo "${BOLD}Extensions to install: ${#EXTENSIONS[@]}${RESET}"
echo ""

# Step 1: Find VS Code CLI
echo "${BOLD}[Step 1/3]${RESET} Locating VS Code CLI..."
echo ""

CODE_BIN=""

# Function to find VS Code CLI
find_vscode_cli() {
    # Try common locations for VS Code CLI
    local candidates=(
        "code"
        "/vscode/vscode-server/bin/*/bin/remote-cli/code"
        "$HOME/.vscode-server/bin/*/bin/remote-cli/code"
    )

    for candidate in "${candidates[@]}"; do
        if [[ "$candidate" == "code" ]]; then
            if command -v code &>/dev/null; then
                echo "$(command -v code)"
                return 0
            fi
        else
            # Expand glob pattern
            for path in $candidate; do
                if [[ -x "$path" ]]; then
                    echo "$path"
                    return 0
                fi
            done
        fi
    done
    return 1
}

# Wait for VS Code CLI to be available (up to 5 minutes)
print_info "Waiting for VS Code CLI to be available (timeout: 5 minutes)..."
WAIT_COUNT=0
for i in {1..300}; do
    if CODE_BIN=$(find_vscode_cli); then
        print_success "Found VS Code CLI: $CODE_BIN"
        break
    fi
    if (( i % 30 == 0 )); then
        WAIT_COUNT=$((WAIT_COUNT + 1))
        print_info "Still waiting... (${i}s elapsed)"
    fi
    sleep 1
done

# Verify code is available
if [[ -z "$CODE_BIN" ]] || ! command -v "$CODE_BIN" &>/dev/null; then
    print_error "VS Code CLI not available after 5 minutes"
    echo ""
    echo "${BOLD}This suggests:${RESET}"
    echo "  • VS Code remote server hasn't fully initialized"
    echo "  • Extensions will be available after VS Code fully loads"
    echo ""
    echo "${BOLD}What to do:${RESET}"
    echo "  Wait a few minutes and try running this script manually:"
    echo "  ${BLUE}bash ${REPO_ROOT}/.devcontainer/scripts/install-extensions-final.sh${RESET}"
    echo ""
    exit 0  # Don't fail - just skip for now
fi

echo ""

# Step 2: Install extensions
echo "${BOLD}[Step 2/3]${RESET} Installing Extensions..."
echo ""

FAILED=0
INSTALLED=0
SKIPPED=0
TOTAL=${#EXTENSIONS[@]}

for ext in "${EXTENSIONS[@]}"; do
    # Check if extension is already installed (case-insensitive)
    if "$CODE_BIN" --list-extensions 2>/dev/null | grep -qi "^${ext}$"; then
        print_info "$ext (already installed)"
        SKIPPED=$((SKIPPED + 1))
    else
        print_info "Installing $ext..."

        # Run installation and capture output
        if install_output=$("$CODE_BIN" --install-extension "$ext" --force 2>&1); then
            # Verify installation by checking extension list
            if "$CODE_BIN" --list-extensions 2>/dev/null | grep -qi "^${ext}$"; then
                print_success "$ext installed successfully"
                INSTALLED=$((INSTALLED + 1))
            else
                print_error "$ext installation failed (not in extension list)"
                FAILED=$((FAILED + 1))
            fi
        else
            print_error "$ext installation failed"
            echo "  ${YELLOW}Error:${RESET} $install_output"
            FAILED=$((FAILED + 1))
        fi
    fi
done

echo ""

# Step 3: Clean up caches
echo "${BOLD}[Step 3/3]${RESET} Cleaning VS Code Caches..."
echo ""

print_info "Removing cached profile data..."
rm -rf "$HOME/.vscode-server/data/CachedProfilesData" 2>/dev/null || true

print_info "Removing cached extension data..."
rm -rf "$HOME/.vscode-server/data/CachedExtensionsData" 2>/dev/null || true

print_success "Cache cleanup complete"
echo ""

# Summary
print_header "Installation Complete"

echo "Completed: ${BOLD}$(date)${RESET}"
echo ""

echo "${BOLD}Installation Summary:${RESET}"
echo "  Total extensions: ${BOLD}$TOTAL${RESET}"
echo "  ${GREEN}Newly installed: $INSTALLED${RESET}"
echo "  ${BLUE}Already installed: $SKIPPED${RESET}"
if [ $FAILED -gt 0 ]; then
    echo "  ${RED}Failed: $FAILED${RESET}"
fi
echo ""

if [ $FAILED -eq 0 ]; then
    print_success "All extensions are now installed!"
    echo ""
    echo "${BOLD}Next step:${RESET}"
    echo "  Reload VS Code window to activate extensions:"
    echo "    1. Press: ${BOLD}Ctrl+Shift+P${RESET}"
    echo "    2. Type: ${BOLD}Developer: Reload Window${RESET}"
    echo "    3. Press: ${BOLD}Enter${RESET}"
else
    print_warning "$FAILED extension(s) failed to install"
    echo ""
    echo "${BOLD}Troubleshooting:${RESET}"
    echo "  • Check log file: ${BLUE}$LOG_FILE${RESET}"
    echo "  • Try reloading VS Code window"
    echo "  • Run this script again if needed"
fi

echo ""
echo "${BOLD}Log file location:${RESET}"
echo "  ${BLUE}$LOG_FILE${RESET}"
echo ""

exit 0

