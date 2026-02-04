#!/usr/bin/env bash
# Install Wolfram Language Server (lsp-wl) and required paclets

set -euo pipefail

REPO_DIR="/home/vscode/.local/share/lsp-wl"
REPO_URL="https://github.com/kenkangxgwe/lsp-wl.git"
KERNEL="/home/vscode/.local/wolfram/engine/14.3/Executables/WolframKernel"

if [[ -d "$REPO_DIR/.git" ]]; then
    git -C "$REPO_DIR" pull --ff-only
else
    mkdir -p "$(dirname "$REPO_DIR")"
    git clone "$REPO_URL" "$REPO_DIR"
fi

"$KERNEL" -noprompt -run 'PacletInstall["CodeParser"]; PacletInstall["CodeInspector"]; Exit[]'

echo "Wolfram Language Server repo: $REPO_DIR"
echo "Paclets installed: CodeParser, CodeInspector"
