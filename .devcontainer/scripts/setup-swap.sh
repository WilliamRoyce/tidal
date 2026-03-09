#!/usr/bin/env bash
# Create swap file for large Wolfram derivations (e.g. Einstein-Maxwell).
# NOTE: swapon does NOT work inside Docker containers on WSL2 (kernel returns
# EINVAL). For Docker/WSL2, increase Docker Desktop memory instead:
#   Docker Desktop → Settings → Resources → Memory → 16 GB+
# This script is kept for non-Docker environments (native Linux, etc.).
set -e

SWAP_SIZE="${SWAP_SIZE:-4G}"
SWAP_FILE="/swapfile"

if swapon --show | grep -q "$SWAP_FILE"; then
    echo "Swap already active: $(swapon --show)"
    exit 0
fi

if ! sudo swapon --show >/dev/null 2>&1; then
    echo "Warning: swapon not available (missing SYS_ADMIN capability). Skipping swap setup."
    exit 0
fi

echo "Creating ${SWAP_SIZE} swap file..."
sudo fallocate -l "$SWAP_SIZE" "$SWAP_FILE"
sudo chmod 600 "$SWAP_FILE"
sudo mkswap "$SWAP_FILE"
if sudo swapon "$SWAP_FILE" 2>/dev/null; then
    echo "Swap enabled: $(free -h | grep Swap)"
else
    echo "Warning: swapon failed (missing SYS_ADMIN capability). Skipping swap setup."
    sudo rm -f "$SWAP_FILE"
fi
