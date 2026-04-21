#!/usr/bin/env bash
# Install PolyChord (Will Handley's nested sampling package).
# Requires gfortran. On HPC with MPI, omit --no-mpi.
#
# Usage:
#   bash scripts/install_polychord.sh          # devcontainer (no MPI)
#   bash scripts/install_polychord.sh --mpi    # HPC with MPI

set -euo pipefail

CLONE_DIR="/tmp/PolyChordLite"
MPI_FLAG="--no-mpi"

if [[ "${1:-}" == "--mpi" ]]; then
    MPI_FLAG=""
    echo "Building PolyChord WITH MPI support"
else
    echo "Building PolyChord WITHOUT MPI (serial mode)"
fi

# Check gfortran
if ! command -v gfortran &>/dev/null; then
    echo "ERROR: gfortran not found. Install with: sudo apt-get install -y gfortran"
    exit 1
fi

# Clone or update
if [[ -d "$CLONE_DIR" ]]; then
    echo "Updating existing clone at $CLONE_DIR"
    cd "$CLONE_DIR" && git pull --ff-only
else
    echo "Cloning PolyChordLite to $CLONE_DIR"
    git clone --depth 1 https://github.com/PolyChord/PolyChordLite.git "$CLONE_DIR"
fi

cd "$CLONE_DIR"

# Ensure setuptools is available for the build
uv pip install setuptools

# Install with appropriate MPI setting
if [[ -n "$MPI_FLAG" ]]; then
    uv pip install --no-build-isolation --config-settings="--global-option=$MPI_FLAG" .
else
    uv pip install --no-build-isolation .
fi

# Verify installation
echo ""
echo "Verifying installation..."
if uv run python -c 'from pypolychord import run_polychord; print("pypolychord OK")'; then
    echo "PolyChord installed successfully."
else
    echo "ERROR: PolyChord installation failed — import check did not pass."
    exit 1
fi
