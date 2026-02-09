#!/usr/bin/env bash
# notify-install-extensions.sh
# VS Code Extension Installation Notification Creator
# Creates a notification file with instructions for manual extension installation
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

print_info() {
    echo "${BLUE}ℹ${RESET} $1"
}

# Configuration
NOTIFICATION_FILE="/home/vscode/.vscode-server/data/extension-install-needed.txt"

echo ""
echo "${BOLD}${BLUE}VS Code Extension Installation Notification${RESET}"
echo "=============================================="
echo ""
echo "This script creates a notification file for the user to manually"
echo "install VS Code extensions after the container is fully initialized."
echo ""

# Create notification file
print_info "Creating notification file..."
mkdir -p "$(dirname "$NOTIFICATION_FILE")"

cat > "$NOTIFICATION_FILE" << 'EOF'
================================================================================
                    VS CODE EXTENSION INSTALLATION REQUIRED
================================================================================

VS Code extensions need to be installed manually due to timing constraints
during container initialization.

TO INSTALL ALL REQUIRED EXTENSIONS:

  1. Open a terminal in VS Code
  2. Run the installation script:
     bash /workspaces/torsion-gertsenshtein/.devcontainer/scripts/install-extensions-final.sh
  3. After completion, reload the VS Code window:
     Press: Ctrl+Shift+P
     Type: "Developer: Reload Window"
     Press: Enter

This only needs to be done once per container rebuild.

EXTENSIONS TO BE INSTALLED:

  Development Tools:
    • ms-python.python          - Python language support
    • ms-python.debugpy          - Python debugging
    • charliermarsh.ruff         - Fast Python linter
    • esbenp.prettier-vscode     - Code formatter

  Utilities:
    • streetsidesoftware.code-spell-checker - Spell checking
    • tamasfe.even-better-toml   - TOML file support
    • peakchen90.open-html-in-browser - HTML preview

  GitHub Integration:
    • github.vscode-github-actions - GitHub Actions support
    • github.copilot             - GitHub Copilot AI
    • github.copilot-chat        - GitHub Copilot Chat
    • github.vscode-pull-request-github - PR management

  AI Assistant:
    • Anthropic.claude-code      - Claude AI integration

================================================================================
EOF

print_success "Notification file created at: $NOTIFICATION_FILE"
echo ""

print_info "The file will be visible to the user upon container startup"
print_info "User can run the installation script manually when ready"
echo ""

exit 0