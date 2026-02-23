#!/usr/bin/env bash
# Proca with Lorentzian Background — Full derive → inspect → simulate → measure
#
# Physics: Two massive vector fields (A, B) in 2+1D coupled via a Lorentzian
# scalar background G(x,y) = g0 / (1 + r^2/R^2).  The A_0 and B_0 components
# are constraints; A_1, A_2, B_1, B_2 are dynamical.
#
# The position-dependent coupling means the constraint solver must handle
# non-uniform source terms.  The Lorentzian profile has algebraic tails,
# testing coefficient evaluation for slowly-decaying backgrounds.
#
# Running this script:
#   cd examples/proca_background && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/proca_background.json
#   tidal simulate ../data/proca_background.json \
#     --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
#     --param g0=1.0 --param R=8.0 \
#     --ic gaussian --ic-component A_1 --ic-amplitude 0.5 --ic-width 3.0 \
#     --grid-shape 64 --bounds=-30:30,-30:30 --t-end 20.0 \
#     --bc periodic,periodic \
#     --output ../data/proca_background_output
#   tidal plot ../data/proca_background_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system
tidal inspect ../data/proca_background.json

# Step 3: Run simulation
# Gaussian IC in A_1 with periodic BCs; constraint solver auto-detects A_0, B_0
tidal simulate ../data/proca_background.json \
  --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 --param g0=1.0 --param R=8.0 \
  --ic gaussian --ic-component A_1 --ic-amplitude 0.5 --ic-width 3.0 \
  --grid-shape 64 --bounds=-30:30,-30:30 --t-end 20.0 \
  --bc periodic,periodic \
  --output ../data/proca_background_output

# Step 4: Measure conversion between vector field groups
# NOTE: Spectral conversion P(k,t) and dispersion omega(k) are NOT available
# for this system because the position-dependent Lorentzian background breaks
# translation invariance.  Only real-space measurements are valid.
tidal measure ../data/proca_background_output \
  --what conversion \
  --source A_0,A_1,A_2 --target B_0,B_1,B_2

# Step 5: Individual plots (saved into the simulation output directory)
tidal plot ../data/proca_background_output --type amplitude --quiet
tidal plot ../data/proca_background_output --type snapshot --field A_1 --time-index -1 --quiet
