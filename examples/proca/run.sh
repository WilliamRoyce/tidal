#!/usr/bin/env bash
# Proca (massive vector) 1+1D — Full derive → inspect → simulate → plot pipeline
#
# To run manually:  cd examples/proca && tidal derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/proca_1d.json

# Run simulation (Gaussian pulse in A_1, Lorenz gauge automatic for m != 0)
tidal simulate ../data/proca_1d.json \
  --param procaMassSquared=1.0 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component A_1 \
  --ic-width 5.0 \
  --t-end 30.0 \
  --scheme scipy \
  --output ../data/proca_output

# Visualize results
tidal plot ../data/proca_output --type heatmap --field A_1 --output ../data/proca_heatmap.png --quiet
tidal plot ../data/proca_output --type amplitude --output ../data/proca_amplitude.png --quiet
tidal plot ../data/proca_output --type energy --output ../data/proca_energy.png --quiet
tidal plot ../data/proca_output --type profile --field A_1 --output ../data/proca_profile.png --quiet
