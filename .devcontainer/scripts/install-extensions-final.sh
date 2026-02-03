#!/bin/bash
# Install VS Code extensions after VS Code server has fully initialized
# This runs via postAttachCommand after VS Code connects

LOG_FILE="/home/vscode/.vscode-server/data/logs/install-extensions.log"
mkdir -p "$(dirname "$LOG_FILE")"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "================================================"
echo "Installing VS Code Extensions"
echo "Started: $(date)"
echo "================================================"

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
)

# The code CLI should be available now that VS Code is running
CODE_BIN=""

# Function to find VS Code CLI
find_vscode_cli() {
    # Try common locations for VS Code CLI
    local candidates=(
        "code"
        "/vscode/vscode-server/bin/*/bin/remote-cli/code"
        "/home/vscode/.vscode-server/bin/*/bin/remote-cli/code"
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
echo "Waiting for VS Code CLI to be available..."
for i in {1..300}; do
    if CODE_BIN=$(find_vscode_cli); then
        echo "✓ Found VS Code CLI: $CODE_BIN"
        break
    fi
    if (( i % 30 == 0 )); then
        echo "  Still waiting... (${i}s elapsed)"
    fi
    sleep 1
done

# Verify code is available
if [[ -z "$CODE_BIN" ]] || ! command -v "$CODE_BIN" &>/dev/null; then
    echo "❌ VS Code CLI not available after 5 minutes"
    echo "This suggests VS Code remote server hasn't fully initialized."
    echo "Extensions may be available after VS Code fully loads."
    echo "Try running this script manually later:"
    echo "  bash /workspaces/torsion-gertsenshtein/.devcontainer/scripts/install-extensions-final.sh"
    exit 0  # Don't fail - just skip for now
fi

echo "Using VS Code CLI: $CODE_BIN"
echo ""

# Install extensions
FAILED=0
for ext in "${EXTENSIONS[@]}"; do
    # Check if extension is installed (case-insensitive)
    if "$CODE_BIN" --list-extensions 2>/dev/null | grep -qi "^${ext}$"; then
        echo "✓ $ext (already installed)"
    else
        echo "Installing $ext..."
        # Run installation and capture full output
        install_output=$("$CODE_BIN" --install-extension "$ext" --force 2>&1)
        echo "$install_output"
        
        # Check if installation was successful by checking if extension now appears in list
        if "$CODE_BIN" --list-extensions 2>/dev/null | grep -qi "^${ext}$"; then
            echo "  ✅ Installation verified - $ext is now installed"
        else
            echo "  ❌ Installation failed - $ext not found in extension list"
            FAILED=$((FAILED + 1))
        fi
    fi
done

echo ""
echo "================================================"
echo "Clearing VS Code caches..."
rm -rf /home/vscode/.vscode-server/data/CachedProfilesData 2>/dev/null || true
rm -rf /home/vscode/.vscode-server/data/CachedExtensionsData 2>/dev/null || true

echo ""
echo "Completed: $(date)"
echo "================================================"
echo ""
echo "✅ Extension installation complete!"
echo ""
echo "If extensions don't appear in the sidebar:"
echo "  Press: Ctrl+Shift+P"
echo "  Type: Developer: Reload Window"
echo "  Press: Enter"
echo ""
echo "================================================"

exit 0

