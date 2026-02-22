#!/usr/bin/env bash
# Klein-Gordon 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: Massive Klein-Gordon scalar field in 3+1D flat spacetime.
# ∂²φ/∂t² = ∇²φ - m²φ with 3D isotropic Gaussian initial condition.
#
# NOTE: 3D data — tidal plot is for 1D/2D only.
#
# Running this script:
#   cd examples/scalar_field_3d && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/klein_gordon_3d.json
#   tidal simulate ../data/klein_gordon_3d.json --param m2=1.0 \
#     --grid-shape 32 --bounds 0:20 --periodic --ic gaussian \
#     --ic-width 2.0 --t-end 8.0

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
  --t-end 8.0
