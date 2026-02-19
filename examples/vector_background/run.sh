#!/usr/bin/env bash
# Vector Background Domain Wall — Full derive → inspect → simulate → measure
#
# Physics: Scalar phi + vector A coupled via an external vector field B_a(x,y)
# where B = (0, 0, B0*tanh(x/W)*exp(-(x^2+y^2)/(2R^2))).  The coupling is a
# localized tanh: antisymmetric across x=0, decaying to ~0 at periodic edges.
#
# This is the first test of the vector background ComponentValue mechanism.
# The generated WLS includes ValidateNoUnresolvedBackgrounds to catch any
# silent failure in the ComponentValue + ToBasis substitution.
#
# NOTE: Spectral conversion P(k,t) and dispersion omega(k) are NOT available
# for this system because the position-dependent background breaks translation
# invariance.  Only real-space measurements (P(t), energy) are valid.
#
# Running this script:
#   cd examples/vector_background && uv run bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system
tidal inspect ../data/vector_background.json

# Step 3: Run simulation
# Gaussian pulse in phi at x=-15 (left of domain wall); disperses into coupling region
tidal simulate ../data/vector_background.json \
  --param mPhi2=1.0 --param mA2=2.0 --param gBV=0.5 \
  --param B0=1.0 --param W=3.0 --param R=8.0 \
  --ic gaussian --ic-width 3.0 --ic-center -15.0,0.0 \
  --grid-shape 64 --bounds -30:30,-30:30 --t-end 20.0 \
  --bc periodic,periodic --scheme scipy \
  --output ../data/vector_background_output

# Step 4: Measure conversion (real-space only — no spectral)
tidal measure ../data/vector_background_output \
  --what conversion \
  --source phi_0 --target A_0,A_1,A_2

# Step 5: Individual plots (saved into the simulation output directory)
tidal plot ../data/vector_background_output --type amplitude --quiet
tidal plot ../data/vector_background_output --type snapshot --time-index -1 --quiet
