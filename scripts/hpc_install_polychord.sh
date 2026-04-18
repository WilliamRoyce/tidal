#!/usr/bin/env bash
# hpc_install_polychord.sh — one-time PolyChord install on CSD3.
#
# Run AFTER `ssh csd3` has a live ControlMaster (authenticated session).
# Executes the install remotely, no interactive prompts needed.
#
# Usage:
#   bash scripts/hpc_install_polychord.sh

set -euo pipefail

readonly HOST="${HPC_HOST:-csd3}"
readonly REMOTE_ROOT="${HPC_ROOT:-/rds/user/wr286/hpc-work/tidal}"

# Verify ControlMaster is live (don't trigger fresh MFA prompt)
if ! ssh -O check "$HOST" 2>/dev/null; then
    cat >&2 <<EOF
error: no live SSH ControlMaster for '$HOST'.

Open one interactively first (password + TOTP):

    ssh $HOST

Then re-run this command.
EOF
    exit 2
fi

echo "==> installing PolyChord on ${HOST}:${REMOTE_ROOT}"

# Remote script: clone, build with MPI, install into the project venv
ssh "$HOST" 'bash -s' <<REMOTE_EOF
set -euo pipefail

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl
module load python/3.11.0-icl

REMOTE_ROOT="${REMOTE_ROOT}"
POLYCHORD_DIR="\$REMOTE_ROOT/.polychord_build"

# Clone or update
if [[ -d "\$POLYCHORD_DIR" ]]; then
    echo "updating existing clone at \$POLYCHORD_DIR"
    cd "\$POLYCHORD_DIR" && git pull --ff-only || true
else
    echo "cloning PolyChordLite to \$POLYCHORD_DIR"
    git clone --depth 1 https://github.com/PolyChord/PolyChordLite.git "\$POLYCHORD_DIR"
fi

cd "\$POLYCHORD_DIR"

# Activate project venv and ensure setuptools is present
source "\$REMOTE_ROOT/.venv/bin/activate"
pip install setuptools

# Build WITH MPI on HPC (serial mode for devcontainer uses --no-mpi)
pip install --no-build-isolation .

# Verify
python -c 'from pypolychord import run_polychord; print("pypolychord installed OK")'
python -c 'from anesthetic import read_chains; print("anesthetic installed OK")'
REMOTE_EOF

echo "==> PolyChord installed on HPC"
