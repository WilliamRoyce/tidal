#!/bin/bash
# Save current Wolfram activation to persistent storage

echo "================================================"
echo "🔄 Wolfram Activation Backup Manager"
echo "================================================"
echo ""

USERBASE_DIR="/home/vscode/.local/wolfram/userbase"
BACKUP_DIR="$USERBASE_DIR/.activation_backup"
CACHE_DIR="/home/vscode/.cache/Wolfram"

if [[ ! -d "$USERBASE_DIR" ]]; then
    echo "❌ ERROR: Wolfram userbase not found at $USERBASE_DIR"
    echo "   Make sure Wolfram Engine is properly mounted"
    exit 1
fi

# Function to backup activation
backup_activation() {
    if [[ -d "$CACHE_DIR/WolframScript" ]]; then
        mkdir -p "$BACKUP_DIR"
        cp -r "$CACHE_DIR"/* "$BACKUP_DIR/"
        echo "✅ Activation backed up to: $BACKUP_DIR"
        echo "   Backup contains:"
        find "$BACKUP_DIR" -type f | sed 's|^|     |'
        echo ""
        echo "   This backup will survive container rebuilds!"
        return 0
    else
        echo "❌ No activation data found in $CACHE_DIR"
        return 1
    fi
}

# Function to restore activation
restore_activation() {
    if [[ -d "$BACKUP_DIR" ]]; then
        mkdir -p "$CACHE_DIR"
        cp -r "$BACKUP_DIR"/* "$CACHE_DIR/"
        echo "✅ Activation restored from backup"
        echo "   Restored files:"
        find "$CACHE_DIR" -type f | sed 's|^|     |'
        return 0
    else
        echo "⚠️  No activation backup found at $BACKUP_DIR"
        return 1
    fi
}

# Function to check status
check_status() {
    echo "📋 Current Status:"
    echo "─────────────────"
    
    echo -n "  Cache directory mounted: "
    if mountpoint -q "$CACHE_DIR" 2>/dev/null; then
        echo "✅ YES (will persist through rebuilds)"
    else
        echo "❌ NO (rebuild required for persistence)"
    fi
    
    echo -n "  Activation data in cache: "
    if [[ -d "$CACHE_DIR/WolframScript" ]] && [[ -n "$(find "$CACHE_DIR/WolframScript" -name "connection_*" -type f 2>/dev/null)" ]]; then
        echo "✅ YES"
    else
        echo "❌ NO"
    fi
    
    echo -n "  Backup exists: "
    if [[ -d "$BACKUP_DIR/WolframScript" ]]; then
        echo "✅ YES"
        echo "     Location: $BACKUP_DIR"
        echo "     Files: $(find "$BACKUP_DIR" -type f | wc -l) files"
    else
        echo "❌ NO"
    fi
    
    echo -n "  WolframScript works: "
    if timeout 10s wolframscript -code "2+2" 2>/dev/null | grep -q "4"; then
        echo "✅ YES"
    else
        echo "❌ NO (needs activation)"
    fi
}

case "${1:-status}" in
    "backup")
        echo "🔄 Backing up current activation..."
        echo ""
        backup_activation
        ;;
    
    "restore")
        echo "🔄 Restoring activation from backup..."
        echo ""
        restore_activation
        ;;
    
    "status"|"check"|"")
        check_status
        ;;
    
    "help"|"-h"|"--help")
        echo "Usage: $0 [backup|restore|status]"
        echo ""
        echo "Commands:"
        echo "  backup   - Save current activation to persistent storage"
        echo "  restore  - Restore activation from persistent storage"
        echo "  status   - Check current activation status (default)"
        echo "  help     - Show this help"
        echo ""
        echo "Background:"
        echo "  Wolfram activation is stored in ~/.cache/Wolfram which is lost on rebuild"
        echo "  unless properly mounted. This script backs up activation data to the"
        echo "  mounted userbase directory so it persists across rebuilds."
        ;;
    
    *)
        echo "❌ Unknown command: $1"
        echo "   Use '$0 help' for usage information"
        exit 1
        ;;
esac

echo ""
echo "================================================"

if [[ "${1:-status}" == "status" ]]; then
    echo ""
    echo "💡 Tips:"
    echo "  • If cache not mounted, rebuild container to apply cache mount"
    echo "  • If no activation, run: wolframscript -activate"
    echo "  • Backup activation after successful activation: $0 backup"
    echo "  • Check diagnostics: bash .devcontainer/check-wolfram.sh"
    echo "  • Full guide: cat .devcontainer/docs/WOLFRAM_GUIDE.md"
fi