#!/bin/bash
# xPerm MathLink Compiler and Installer
# Compiles xPerm from source with current GLIBC and installs with proper library path

set -e

XPERM_DIR="$WOLFRAM_USERBASE/Applications/xAct/xPerm/mathlink"
MATHLINK_DIR="/home/vscode/.local/wolfram/engine/14.3/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions"

echo "=== xPerm MathLink Compiler and Installer ==="
echo ""

# Check prerequisites
echo "Checking prerequisites..."
if ! command -v gcc &> /dev/null; then
    echo "❌ gcc not found. Installing build tools..."
    sudo apt update && sudo apt install -y build-essential uuid-dev
fi

if [[ ! -f "$MATHLINK_DIR/mprep" ]]; then
    echo "❌ Wolfram MathLink SDK not found at $MATHLINK_DIR"
    exit 1
fi

cd "$XPERM_DIR"

echo "✅ Prerequisites satisfied"
echo ""

# Backup original files
echo "Creating backups..."
if [[ ! -f "xperm.c.original" ]]; then
    cp xperm.c xperm.c.original
fi

if [[ -f "xperm.linux.64-bit" && ! -f "xperm.linux.64-bit.factory" ]]; then
    cp xperm.linux.64-bit xperm.linux.64-bit.factory
fi

echo "✅ Backups created"
echo ""

# Generate MathLink template
echo "Generating MathLink template..."
"$MATHLINK_DIR/mprep" xperm.tm -o xpermp.c
echo "✅ Template generated"

# Compile xPerm
echo "Compiling xPerm with current GLIBC..."
gcc -I"$MATHLINK_DIR" -L"$MATHLINK_DIR" -O2 xpermp.c \
    -lML64i4 -lpthread -lrt -lstdc++ -ldl -luuid -lm \
    -o xperm.linux.64-bit.compiled

echo "✅ Compilation successful"

# Create wrapper script
echo "Creating MathLink wrapper..."
cat > xperm.linux.64-bit.wrapper << 'EOF'
#!/bin/bash
# xPerm MathLink Wrapper - Auto-generated
MATHLINK_DIR="/home/vscode/.local/wolfram/engine/14.3/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions"
export LD_LIBRARY_PATH="$MATHLINK_DIR:$LD_LIBRARY_PATH"
exec "$(dirname "$0")/xperm.linux.64-bit.compiled" "$@"
EOF

chmod +x xperm.linux.64-bit.wrapper
echo "✅ Wrapper created"

# Install new version
echo "Installing new xPerm binary..."
if [[ -f "xperm.linux.64-bit" ]]; then
    mv xperm.linux.64-bit xperm.linux.64-bit.old
fi
mv xperm.linux.64-bit.wrapper xperm.linux.64-bit

echo "✅ Installation complete"
echo ""

# Test installation
echo "Testing xPerm MathLink connection..."
if timeout 10 wolframscript -code '<<xAct`xPerm`; Print["Connection test: ", $xpermQ];' 2>/dev/null | grep -q "Connection established"; then
    echo "✅ xPerm MathLink working perfectly!"
    echo ""
    echo "🎉 SUCCESS: xPerm advanced algorithms now available!"
    echo "   - Fast permutation group computations"
    echo "   - Strong generating sets"
    echo "   - Stabilizer chain algorithms"
    echo "   - All xPerm functions at full performance"
else
    echo "⚠️  Installation completed but test failed. Manual verification needed."
fi

echo ""
echo "Files created/modified:"
echo "  - xperm.linux.64-bit (wrapper script)"
echo "  - xperm.linux.64-bit.compiled (new binary)"
echo "  - xperm.c.original (backup)"
echo "  - xperm.linux.64-bit.factory (factory backup)"
echo ""
echo "To restore factory version: mv xperm.linux.64-bit.factory xperm.linux.64-bit"