#!/bin/bash
# Quick Wolfram Engine health check

echo "🔬 Wolfram Engine Health Check"
echo "================================"

# Check wolframscript
if ! command -v wolframscript &>/dev/null; then
    echo "❌ wolframscript not found"
    exit 1
fi

echo "✓ wolframscript: $(command -v wolframscript)"

# Test WolframKernel
KERNEL_TEST=$(/home/vscode/.local/wolfram/engine/14.3/Executables/WolframKernel -noprompt -run "Print[2+2]; Exit[]" 2>&1)
if [[ "$KERNEL_TEST" == *"4"* ]]; then
    echo "✅ WolframKernel: Working"
else
    echo "❌ WolframKernel: Failed"
    echo "   Error: $KERNEL_TEST"
fi

# Test WolframScript
SCRIPT_TEST=$(wolframscript -code "3+3" 2>&1)
if [[ "$SCRIPT_TEST" == "6" ]]; then
    echo "✅ WolframScript: Working"
else
    echo "❌ WolframScript: Failed"
    echo "   Error: $SCRIPT_TEST"
    echo "   Try: bash .devcontainer/wolfram-activation-manager.sh restore"
fi

# Check persistence
echo ""
echo "📁 Persistence Status:"
if mountpoint -q /home/vscode/.cache/Wolfram 2>/dev/null; then
    echo "✅ Cache mounted (activation will persist)"
else
    echo "⚠️  Cache not mounted (rebuild required)"
fi

if [[ -d "/home/vscode/.local/wolfram/userbase/.activation_backup" ]] && [[ -n "$(find "/home/vscode/.local/wolfram/userbase/.activation_backup" -name "connection_*" -type f 2>/dev/null)" ]]; then
    echo "✅ Activation backup exists"
else
    echo "⚠️  No activation backup found"
fi

echo ""
echo "📚 Help:"
echo "  Full guide: cat .devcontainer/WOLFRAM_GUIDE.md"
echo "  Manage activation: bash .devcontainer/wolfram-activation-manager.sh"