#!/bin/bash
# xPerm MathLink Compatibility Fix
# Patches xPerm to work without requiring GLIBC 2.38

XPERM_FILE="$WOLFRAM_USERBASE/Applications/xAct/xPerm/xPerm.m"

if [[ -f "$XPERM_FILE" ]]; then
    echo "Applying xPerm compatibility fix..."
    
    # Create backup if it doesn't exist
    if [[ ! -f "$XPERM_FILE.backup" ]]; then
        cp "$XPERM_FILE" "$XPERM_FILE.backup"
        echo "Created backup: xPerm.m.backup"
    fi
    
    # Replace the problematic mathlink connection logic
    sed -i '/Print\[If\[result,"Connection established\.","Connection failed\."\]\]/c\
    If[FileExistsQ[$xpermExecutable],\
        Check[$xpermLink=Install[$xpermExecutable],result=False];\
        Print[If[result,"Connection established.","Connection failed."]],\
        Print["MathLink binary not compatible. Using pure Mathematica mode."];\
        result=False\
    ]' "$XPERM_FILE"
    
    echo "✅ xPerm patched for compatibility"
else
    echo "❌ xPerm.m not found at $XPERM_FILE"
fi