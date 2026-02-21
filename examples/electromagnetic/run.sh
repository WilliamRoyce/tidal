#!/usr/bin/env bash
# Electromagnetic 2+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: Maxwell's equations in 2+1D. A_0 is a constraint (Gauss's law,
# time_order=0), A_1 and A_2 are dynamical. For a plane wave propagating
# along x, the transverse polarization is A_2 (= A_y). The constraint solver
# ensures Gauss's law (∇²A_0 = ∂_x π_1 + ∂_y π_2) is satisfied at each step.
#
# NOTE: 1+1D EM has zero propagating DOF (no transverse dimension), so
# the example requires 2+1D. Exciting A_2 with --ic plane-wave gives a
# physical EM wave with E_y and B_z oscillating at c = 1.
#
# Running this script:
#   cd examples/electromagnetic && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/em_3d.json
#   tidal simulate ../data/em_3d.json --ic plane-wave --ic-component A_2 \
#     --grid-shape 64 --bounds 0:50 --periodic --t-end 25.0 --scheme scipy \
#     --output ../data/em_output
#   tidal plot ../data/em_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/em_3d.json

# Run simulation (plane wave in A_2 — transverse polarization for x-propagation)
# A_0 is a constraint (time_order=0), A_1 and A_2 are dynamical
tidal simulate ../data/em_3d.json \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic plane-wave \
  --ic-component A_2 \
  --t-end 25.0 \
  --scheme scipy \
  --output ../data/em_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/em_output --type heatmap --field A_2 --quiet
tidal plot ../data/em_output --type amplitude --quiet
tidal plot ../data/em_output --type snapshot --time-index -1 --quiet
