#!/usr/bin/env bash
# Electromagnetic 1+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: Maxwell's equations in 1+1D. A_0 is a constraint (Gauss's law,
# time_order=0), only A_1 is dynamical. Plane-wave IC provides non-zero
# conjugate momentum π_1 (= electric field), giving massless propagation at c=1.
#
# NOTE: Gaussian IC with π=0 is a pure-gauge configuration (zero E-field)
# and will NOT propagate. Always use --ic plane-wave for EM.
#
# Running this script:
#   cd examples/electromagnetic && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/em_1d.json
#   tidal simulate ../data/em_1d.json --ic plane-wave --ic-component A_1 \
#     --grid-shape 256 --bounds 0:100 --periodic --t-end 25.0 --scheme scipy \
#     --output ../data/em_output
#   tidal plot ../data/em_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/em_1d.json

# Run simulation (plane wave in A_1 component — provides non-zero E-field)
# A_0 is a constraint (time_order=0), only A_1 is dynamical
tidal simulate ../data/em_1d.json \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-component A_1 \
  --t-end 25.0 \
  --scheme scipy \
  --output ../data/em_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/em_output --type heatmap --field A_1 --quiet
tidal plot ../data/em_output --type amplitude --quiet
tidal plot ../data/em_output --type snapshot --time-index -1 --quiet
