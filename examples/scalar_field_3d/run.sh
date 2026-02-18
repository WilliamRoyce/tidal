#!/usr/bin/env bash
# CLI equivalents for the Klein-Gordon 3+1D example
# See also: kg_3d_simulation.py (Python simulation)
#
# To run manually:  cd examples/scalar_field_3d && tidal derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/klein_gordon_3d.json

# Run simulation (3D Gaussian pulse at domain center)
tidal simulate ../data/klein_gordon_3d.json \
  --param m2=1.0 \
  --grid-shape 32 \
  --bounds 0:20 \
  --periodic \
  --ic gaussian \
  --ic-width 2.0 \
  --t-end 8.0 \
  --scheme scipy
