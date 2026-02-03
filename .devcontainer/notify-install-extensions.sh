#!/bin/bash
# Auto-installer trigger that creates a notification file
# This runs immediately and creates instructions for manual installation

NOTIFICATION_FILE="/home/vscode/.vscode-server/data/extension-install-needed.txt"

cat > "$NOTIFICATION_FILE" << 'EOF'
========================================
🔧 EXTENSION INSTALLATION REQUIRED
========================================

VS Code extensions need to be installed manually due to timing constraints.

TO INSTALL ALL REQUIRED EXTENSIONS:
1. Open a terminal in VS Code
2. Run: bash /workspaces/torsion-gertsenshtein/.devcontainer/install-extensions-final.sh
3. After completion, reload window: Ctrl+Shift+P → "Developer: Reload Window"

This only needs to be done once per rebuild.

Extensions to be installed:
- Python support (ms-python.python + dependencies)
- Ruff linting (charliermarsh.ruff)
- Prettier formatting (esbenp.prettier-vscode)
- Spell checking (streetsidesoftware.code-spell-checker)
- TOML support (tamasfe.even-better-toml)
- HTML preview (peakchen90.open-html-in-browser)
- GitHub Actions (github.vscode-github-actions)
- GitHub Copilot (github.copilot + copilot-chat)
- GitHub PR management (github.vscode-pull-request-github)
- Claude AI (anthropic.claude-code)

========================================
EOF

echo "Extension installation notification created at: $NOTIFICATION_FILE"
echo "User will be prompted to run installation script manually."

exit 0