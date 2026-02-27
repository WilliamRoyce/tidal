#!/usr/bin/env bash
# reindex-claude-sessions.sh -- Rebuild sessions-index.json and history.jsonl
# from orphaned .jsonl session files that the VS Code extension fails to index.
# Workaround for: https://github.com/anthropics/claude-code/issues/18619
#
# Usage: bash .devcontainer/scripts/reindex-claude-sessions.sh
set -euo pipefail

CLAUDE_DIR="$HOME/.claude"
PROJECT_DIR="$CLAUDE_DIR/projects/-workspaces-torsion-gertsenshtein"
INDEX_FILE="$PROJECT_DIR/sessions-index.json"
HISTORY_FILE="$CLAUDE_DIR/history.jsonl"

if [ ! -d "$PROJECT_DIR" ]; then
  echo "No project directory found: $PROJECT_DIR"
  exit 0
fi

# Count session files (exclude directories like memory/)
session_count=$(find "$PROJECT_DIR" -maxdepth 1 -name "*.jsonl" -type f 2>/dev/null | wc -l)
if [ "$session_count" -eq 0 ]; then
  echo "No session files found in $PROJECT_DIR"
  exit 0
fi

echo "Found $session_count session file(s). Rebuilding index..."

python3 << 'PYEOF'
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

project_dir = Path(os.path.expanduser("~/.claude/projects/-workspaces-torsion-gertsenshtein"))
index_file = project_dir / "sessions-index.json"
history_file = Path(os.path.expanduser("~/.claude/history.jsonl"))
project_path = "/workspaces/torsion-gertsenshtein"

# Load existing history session IDs
existing_history = set()
if history_file.exists():
    for line in history_file.read_text().strip().split("\n"):
        if not line.strip():
            continue
        try:
            existing_history.add(json.loads(line)["sessionId"])
        except (json.JSONDecodeError, KeyError):
            pass

entries = []
new_history = []

for jsonl_file in sorted(project_dir.glob("*.jsonl")):
    session_id = jsonl_file.stem
    stat = jsonl_file.stat()
    mtime_ms = int(stat.st_mtime * 1000)

    # Parse the session file to extract metadata
    first_prompt = "Untitled"
    message_count = 0
    git_branch = ""
    created_ts = None
    is_sidechain = False

    try:
        with open(jsonl_file) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    d = json.loads(line)
                except json.JSONDecodeError:
                    continue

                msg_type = d.get("type", "")

                if msg_type in ("user", "assistant"):
                    message_count += 1

                # Extract first user prompt
                if msg_type == "user" and first_prompt == "Untitled":
                    msg = d.get("message", {})
                    content = msg.get("content", "") if isinstance(msg, dict) else ""
                    if isinstance(content, list):
                        for block in content:
                            if isinstance(block, dict) and block.get("type") == "text":
                                text = block.get("text", "").strip()
                                if text:
                                    first_prompt = text[:120]
                                    break
                    elif isinstance(content, str) and content.strip():
                        first_prompt = content.strip()[:120]

                # Extract git branch from first message with it
                if not git_branch and d.get("gitBranch"):
                    git_branch = d["gitBranch"]

                # Extract creation timestamp from first entry
                if created_ts is None and d.get("timestamp"):
                    ts = d["timestamp"]
                    if isinstance(ts, str):
                        created_ts = ts
                    elif isinstance(ts, (int, float)):
                        created_ts = datetime.fromtimestamp(
                            ts / 1000 if ts > 1e12 else ts, tz=timezone.utc
                        ).isoformat()

                # Check sidechain
                if d.get("isSidechain"):
                    is_sidechain = True

    except Exception as e:
        print(f"  Warning: error reading {jsonl_file.name}: {e}", file=sys.stderr)
        continue

    if created_ts is None:
        created_ts = datetime.fromtimestamp(stat.st_ctime, tz=timezone.utc).isoformat()

    modified_ts = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()

    entry = {
        "sessionId": session_id,
        "fullPath": str(jsonl_file),
        "fileMtime": mtime_ms,
        "firstPrompt": first_prompt,
        "messageCount": message_count,
        "created": created_ts,
        "modified": modified_ts,
        "gitBranch": git_branch,
        "projectPath": project_path,
        "isSidechain": is_sidechain,
    }
    entries.append(entry)
    print(f"  {session_id[:8]}... | {message_count:3d} msgs | {git_branch or 'no-branch'} | {first_prompt[:60]}")

    # Add to history if missing
    if session_id not in existing_history:
        history_entry = {
            "display": first_prompt,
            "pastedContents": {},
            "timestamp": mtime_ms,
            "project": str(project_dir),
            "sessionId": session_id,
        }
        new_history.append(history_entry)

# Write sessions-index.json
index_data = {"version": 1, "entries": entries}
with open(index_file, "w") as f:
    json.dump(index_data, f, indent=2)
print(f"\nWrote {len(entries)} entries to sessions-index.json")

# Append to history.jsonl
if new_history:
    with open(history_file, "a") as f:
        for entry in new_history:
            f.write(json.dumps(entry) + "\n")
    print(f"Added {len(new_history)} new entries to history.jsonl")
else:
    print("history.jsonl already up to date")

PYEOF

echo "Done. Restart VS Code or reload the window to see sessions in the dropdown."
