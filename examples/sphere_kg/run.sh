#!/usr/bin/env bash
# Sphere Klein-Gordon 2+1D (stereographic projection) — Full pipeline
#
# Physics: Klein-Gordon on S² via stereographic projection. Position-dependent
# conformal factor Ω(r) = 2R²/(R² + r²) produces variable wave speed.
# The metric is coordinate-dependent, exercising Christoffel auto-detection.
#
# Running this script:
#   cd examples/sphere_kg && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/sphere_kg.json
#   tidal simulate ../data/sphere_kg.json --param sphR=2.0 --param sphm2=0.0 \
#     --grid-shape 128 --bounds=-8:8 --periodic --ic gaussian --ic-width 0.8 \
#     --t-end 10.0 --scheme scipy --output ../data/sphere_kg_output
#   tidal plot ../data/sphere_kg_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/sphere_kg.json

# Run simulation (2D Gaussian on stereographic plane)
# Position-dependent wave speed from conformal factor Omega(r) = 2R^2/(R^2 + r^2)
tidal simulate ../data/sphere_kg.json \
  --param sphR=2.0 --param sphm2=0.0 \
  --grid-shape 128 \
  --bounds -8:8 \
  --periodic \
  --ic gaussian \
  --ic-width 0.8 \
  --t-end 10.0 \
  --scheme scipy \
  --output ../data/sphere_kg_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/sphere_kg_output --type snapshot --time-index 0 --quiet
tidal plot ../data/sphere_kg_output --type snapshot --time-index -1 --quiet
tidal plot ../data/sphere_kg_output --type profile --cross-section y=0.0 --quiet
tidal plot ../data/sphere_kg_output --type amplitude --quiet
