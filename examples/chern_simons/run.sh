#!/usr/bin/env bash
# Chern-Simons 2+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: Topological Chern-Simons gauge theory in 2+1D with coupling κ.
# The epsilon tensor (ε^abc) produces parity-violating dynamics. A_0 is a
# constraint (time_order=0), A_1 and A_2 are dynamical.
#
# Running this script:
#   cd examples/chern_simons && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/chern_simons_3d.json
#   tidal simulate ../data/chern_simons_3d.json --param kappa=0.5 \
#     --grid-shape 64 --bounds 0:50 --periodic --ic gaussian \
#     --ic-component A_1 --ic-width 5.0 --t-end 10.0 --scheme scipy \
#     --output ../data/chern_simons_output
#   tidal plot ../data/chern_simons_output --type amplitude --quiet

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
tidal plot ../data/chern_simons_output --type snapshot --field A_0 --time-index 0 --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_1 --time-index 0 --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_2 --time-index 0 --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_0 --time-index -1 --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_1 --time-index -1 --quiet
tidal plot ../data/chern_simons_output --type snapshot --field A_2 --time-index -1 --quiet
tidal plot ../data/chern_simons_output --type amplitude --quiet
