#!/usr/bin/env bash
# sync-claude-memory.sh -- Backup and restore Claude Code memory files
# Usage:
#   sync-claude-memory.sh backup   -- Copy memory -> .claude-memory-backup/
#   sync-claude-memory.sh restore  -- Copy .claude-memory-backup/ -> memory
#   sync-claude-memory.sh status   -- Show sync status
set -euo pipefail

MEMORY_DIR="/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory"
BACKUP_DIR="/workspaces/torsion-gertsenshtein/.claude-memory-backup"

case "${1:-status}" in
  backup)
    if [ ! -d "$MEMORY_DIR" ] || [ -z "$(ls -A "$MEMORY_DIR" 2>/dev/null)" ]; then
      echo "ERROR: Memory directory is empty or missing: $MEMORY_DIR"
      exit 1
    fi
    mkdir -p "$BACKUP_DIR"
    rsync -av --delete "$MEMORY_DIR/" "$BACKUP_DIR/" --exclude='.last-sync'
    date -Iseconds > "$BACKUP_DIR/.last-sync"
    echo "Backup complete: $MEMORY_DIR -> $BACKUP_DIR"
    ;;
  restore)
    if [ ! -d "$BACKUP_DIR" ] || [ -z "$(ls -A "$BACKUP_DIR" 2>/dev/null)" ]; then
      echo "ERROR: Backup directory is empty or missing: $BACKUP_DIR"
      exit 1
    fi
    mkdir -p "$MEMORY_DIR"
    rsync -av "$BACKUP_DIR/" "$MEMORY_DIR/" --exclude='.last-sync'
    echo "Restore complete: $BACKUP_DIR -> $MEMORY_DIR"
    ;;
  status)
    echo "=== Claude Memory Sync Status ==="
    echo ""
    echo "Memory dir: $MEMORY_DIR"
    if [ -d "$MEMORY_DIR" ]; then
      count=$(find "$MEMORY_DIR" -type f | wc -l)
      echo "  Files: $count"
      if [ "$count" -gt 0 ]; then
        find "$MEMORY_DIR" -type f -printf "  - %f (%s bytes)\n"
      fi
    else
      echo "  MISSING"
    fi
    echo ""
    echo "Backup dir: $BACKUP_DIR"
    if [ -d "$BACKUP_DIR" ]; then
      count=$(find "$BACKUP_DIR" -type f -not -name '.last-sync' | wc -l)
      echo "  Files: $count"
      if [ "$count" -gt 0 ]; then
        find "$BACKUP_DIR" -type f -not -name '.last-sync' -printf "  - %f (%s bytes)\n"
      fi
      if [ -f "$BACKUP_DIR/.last-sync" ]; then
        echo "  Last sync: $(cat "$BACKUP_DIR/.last-sync")"
      fi
    else
      echo "  MISSING"
    fi
    ;;
  *)
    echo "Usage: $0 {backup|restore|status}"
    exit 1
    ;;
esac
