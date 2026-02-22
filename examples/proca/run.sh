#!/usr/bin/env bash
# Proca (massive vector) 1+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: Massive vector field in 1+1D Lorenz gauge. The Proca mass term
# (m² A_a A^a) gives dispersive propagation: ω² = k² + m². A_0 is a constraint,
# only A_1 is dynamical.
#
# Running this script:
#   cd examples/proca && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/proca_1d.json
#   tidal simulate ../data/proca_1d.json --param procaMassSquared=1.0 \
#     --grid-shape 512 --bounds 0:100 --periodic --ic gaussian \
#     --ic-component A_1 --ic-width 5.0 --t-end 30.0 \
#     --output ../data/proca_output
#   tidal plot ../data/proca_output --type amplitude --quiet

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
  --output ../data/proca_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/proca_output --type heatmap --field A_1 --quiet
tidal plot ../data/proca_output --type amplitude --quiet
tidal plot ../data/proca_output --type energy --quiet
tidal plot ../data/proca_output --type profile --field A_1 --quiet
