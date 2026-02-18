#!/usr/bin/env bash
# Klein-Gordon 1+1D — Full derive → inspect → simulate → plot pipeline
#
# To run manually:  cd examples/scalar_field && tidal derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Or preview the generated Wolfram script without running:
# tidal derive theory.toml --dry-run

# Inspect the equation system
tidal inspect ../data/klein_gordon_1d.json

# Run simulation (Gaussian pulse, disk-backed output)
tidal simulate ../data/klein_gordon_1d.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-width 5.0 \
  --t-end 30.0 \
  --scheme scipy \
  --snapshots 0.1 \
  --output ../data/scalar_field_output

# Visualize results
tidal plot ../data/scalar_field_output --type heatmap --output ../data/scalar_field_heatmap.png --quiet
tidal plot ../data/scalar_field_output --type amplitude --output ../data/scalar_field_amplitude.png --quiet
tidal plot ../data/scalar_field_output --type snapshot --time-index 0 --output ../data/scalar_field_t0.png --quiet
tidal plot ../data/scalar_field_output --type snapshot --time-index -1 --output ../data/scalar_field_final.png --quiet
tidal plot ../data/scalar_field_output --type profile --output ../data/scalar_field_profile.png --quiet
