#!/usr/bin/env bash
# Sphere Klein-Gordon 2+1D (stereographic projection) — Full pipeline
#
# NOTE: The derive step uses stereographic projection coordinates with
# position-dependent metric.
# The simulate step works fully via CLI (all periodic BCs, 2D Cartesian grid).
#
# To run manually:  cd examples/sphere_kg

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

# Visualize results
tidal plot ../data/sphere_kg_output --type snapshot --time-index 0 --output ../data/sk_initial.png --quiet
tidal plot ../data/sphere_kg_output --type snapshot --time-index -1 --output ../data/sk_final.png --quiet
tidal plot ../data/sphere_kg_output --type profile --cross-section y=0.0 --output ../data/sk_profile.png --quiet
tidal plot ../data/sphere_kg_output --type amplitude --output ../data/sk_amplitude.png --quiet
