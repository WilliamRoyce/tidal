#!/usr/bin/env bash
# Polar Klein-Gordon 2+1D — Full derive → inspect → simulate → plot pipeline
#
# NOTE: The derive step uses polar coordinates (r, theta) with a coordinate-dependent
# metric.
# The simulation works via CLI with --bc and --ic formula flags.
#
# To run manually:  cd examples/polar_kg

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/polar_kg.json

# Run simulation (Gaussian ring at r=3, Neumann in r, periodic in theta)
# Coordinates: x=r, y=theta
tidal simulate ../data/polar_kg.json \
  --param polm2=0.5 \
  --grid-shape 128 \
  --bounds 0.5:10,0:6.283185 \
  --bc neumann,periodic \
  --ic formula \
  --ic-formula "np.exp(-(x - 3.0)**2 / 0.5)" \
  --t-end 8.0 \
  --scheme scipy \
  --output ../data/polar_kg_output

# Visualize results
tidal plot ../data/polar_kg_output --type snapshot --time-index 0 --output ../data/pk_initial.png --quiet
tidal plot ../data/polar_kg_output --type snapshot --time-index -1 --output ../data/pk_final.png --quiet
tidal plot ../data/polar_kg_output --type amplitude --output ../data/pk_amplitude.png --quiet
