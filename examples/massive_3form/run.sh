#!/usr/bin/env bash
# Massive 3-Form 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: Antisymmetric rank-3 tensor: 64 components reduce to 4 independent
# (C_0..C_3) via epsilon symmetry. Each component satisfies a massive
# Klein-Gordon equation: ∂²C_i/∂t² = ∇²C_i - m²C_i.
#
# NOTE: 3D data — tidal plot is for 1D/2D only.
#
# Running this script:
#   cd examples/massive_3form && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/massive_3form.json
#   tidal simulate ../data/massive_3form.json --param m2=1.0 \
#     --grid-shape 16 --bounds 0:10 --periodic --ic gaussian \
#     --ic-component C_0 --ic-width 1.5 --t-end 5.0 --scheme scipy

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/massive_3form.json

# Run simulation (Gaussian pulse in C_0, other components start at zero)
tidal simulate ../data/massive_3form.json \
  --param m2=1.0 \
  --grid-shape 16 \
  --bounds 0:10 \
  --periodic \
  --ic gaussian \
  --ic-component C_0 \
  --ic-width 1.5 \
  --t-end 5.0 \
  --scheme scipy
