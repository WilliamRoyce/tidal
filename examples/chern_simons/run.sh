#!/usr/bin/env bash
# Chern-Simons 2+1D — Full derive → inspect → simulate → plot pipeline
#
# To run manually:  cd examples/chern_simons && tidal derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/chern_simons_3d.json

# Run simulation (Gaussian pulse in A_1 component)
# A_0 is a constraint (time_order=0), A_1 and A_2 are dynamical
tidal simulate ../data/chern_simons_3d.json \
  --param kappa=0.5 \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian \
  --ic-component A_1 \
  --ic-width 5.0 \
  --t-end 10.0 \
  --scheme scipy \
  --output ../data/chern_simons_output

# Visualize results — initial and final snapshots for each field
tidal plot ../data/chern_simons_output --type snapshot --field A_0 --time-index 0 --output ../data/cs_A0_t0.png --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_1 --time-index 0 --output ../data/cs_A1_t0.png --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_2 --time-index 0 --output ../data/cs_A2_t0.png --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_0 --time-index -1 --output ../data/cs_A0_final.png --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_1 --time-index -1 --output ../data/cs_A1_final.png --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_2 --time-index -1 --output ../data/cs_A2_final.png --quiet
tidal plot ../data/chern_simons_output --type amplitude --output ../data/cs_amplitude.png --quiet
