#!/usr/bin/env bash
# sync-claude-memory.sh -- Backup and restore Claude Code memory + plans files
# Usage:
#   sync-claude-memory.sh backup   -- Copy memory + plans -> .claude-memory-backup/ + .claude-plans-backup/
#   sync-claude-memory.sh restore  -- Copy backups -> memory + plans
#   sync-claude-memory.sh status   -- Show sync status
set -euo pipefail

MEMORY_DIR="/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory"
MEMORY_BACKUP="/workspaces/torsion-gertsenshtein/.claude-memory-backup"

PLANS_DIR="/home/vscode/.claude/plans"
PLANS_BACKUP="/workspaces/torsion-gertsenshtein/.claude-plans-backup"

PROJECT_DIR="/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein"
PROJECT_BACKUP="/workspaces/torsion-gertsenshtein/.claude-project-backup"

_backup_dir() {
  local src="$1" dst="$2" label="$3"
  if [ ! -d "$src" ] || [ -z "$(ls -A "$src" 2>/dev/null)" ]; then
    echo "SKIP: $label source is empty or missing: $src"
    return
  fi
  mkdir -p "$dst"
  rsync -av --delete "$src/" "$dst/" --exclude='.last-sync'
  echo "$label backup complete: $src -> $dst"
}

_restore_dir() {
  local src="$1" dst="$2" label="$3"
  if [ ! -d "$src" ] || [ -z "$(ls -A "$src" 2>/dev/null)" ]; then
    echo "SKIP: $label backup is empty or missing: $src"
    return
  fi
  mkdir -p "$dst"
  rsync -av "$src/" "$dst/" --exclude='.last-sync'
  echo "$label restore complete: $src -> $dst"
}

_status_dir() {
  local dir="$1" label="$2"
  echo "$label: $dir"
  if [ -d "$dir" ]; then
    count=$(find "$dir" -type f -not -name '.last-sync' | wc -l)
    echo "  Files: $count"
    if [ "$count" -gt 0 ]; then
      find "$dir" -type f -not -name '.last-sync' -printf "  - %f (%s bytes)\n" | sort
    fi
  else
    echo "  MISSING"
  fi
}

case "${1:-status}" in
  backup)
    _backup_dir "$MEMORY_DIR" "$MEMORY_BACKUP" "Memory"
    _backup_dir "$PLANS_DIR" "$PLANS_BACKUP" "Plans"
    # Back up project-level settings (CLAUDE.md overrides, settings.json, etc.)
    # but NOT session transcripts (*.jsonl) which are huge
    mkdir -p "$PROJECT_BACKUP"
    find "$PROJECT_DIR" -maxdepth 1 -type f -not -name '*.jsonl' -exec cp {} "$PROJECT_BACKUP/" \;
    echo "Project settings backup complete"
    date -Iseconds > "$MEMORY_BACKUP/.last-sync"
    echo ""
    echo "All backups complete at $(date -Iseconds)"
    ;;
  restore)
    _restore_dir "$MEMORY_BACKUP" "$MEMORY_DIR" "Memory"
    _restore_dir "$PLANS_BACKUP" "$PLANS_DIR" "Plans"
    if [ -d "$PROJECT_BACKUP" ] && [ -n "$(ls -A "$PROJECT_BACKUP" 2>/dev/null)" ]; then
      mkdir -p "$PROJECT_DIR"
      cp "$PROJECT_BACKUP"/* "$PROJECT_DIR/" 2>/dev/null || true
      echo "Project settings restore complete"
    fi
    echo ""
    echo "All restores complete"
    ;;
  status)
    echo "=== Claude Memory & Plans Sync Status ==="
    echo ""
    _status_dir "$MEMORY_DIR" "Memory (live)"
    echo ""
    _status_dir "$MEMORY_BACKUP" "Memory (backup)"
    echo ""
    _status_dir "$PLANS_DIR" "Plans (live)"
    echo ""
    _status_dir "$PLANS_BACKUP" "Plans (backup)"
    echo ""
    _status_dir "$PROJECT_DIR" "Project settings (live)"
    echo ""
    _status_dir "$PROJECT_BACKUP" "Project settings (backup)"
    if [ -f "$MEMORY_BACKUP/.last-sync" ]; then
      echo ""
      echo "Last sync: $(cat "$MEMORY_BACKUP/.last-sync")"
    fi
    ;;
  *)
    echo "Usage: $0 {backup|restore|status}"
    exit 1
    ;;
esac
