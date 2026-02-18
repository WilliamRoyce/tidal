#!/usr/bin/env bash
# Electromagnetic 1+1D — Full derive → inspect → simulate → plot pipeline
#
# To run manually:  cd examples/electromagnetic && tidal derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/em_1d.json

# Run simulation (Gaussian pulse in A_1 component)
# A_0 is a constraint (time_order=0), only A_1 is dynamical
tidal simulate ../data/em_1d.json \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component A_1 \
  --ic-width 5.0 \
  --t-end 25.0 \
  --scheme scipy \
  --output ../data/em_output

# Visualize results
tidal plot ../data/em_output --type heatmap --field A_1 --output ../data/em_heatmap.png --quiet
tidal plot ../data/em_output --type amplitude --output ../data/em_amplitude.png --quiet
tidal plot ../data/em_output --type snapshot --time-index -1 --output ../data/em_final.png --quiet
