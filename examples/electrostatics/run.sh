#!/usr/bin/env bash
# Electrostatics 1+1D — Full derive → inspect → simulate pipeline
#
# Physics: Free scalar field wave equation derived from a Lagrangian,
# demonstrating constraint-mode (Laplace equation) solving. The wave
# equation is ∂²φ/∂t² = ∂²φ/∂x².
#
# Running this script:
#   cd examples/electrostatics && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/electrostatics_laplace.json
#   tidal simulate ../data/electrostatics_laplace.json \
#     --grid-shape 64 --bounds=-5:5 --t-end 5.0 --ic gaussian --no-plot

set -euo pipefail
cd "$(dirname "$0")"

# Derive wave equation from TOML config (Lagrangian-derived)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/electrostatics_laplace.json

# Simulate: wave equation with Gaussian initial condition
tidal simulate ../data/electrostatics_laplace.json \
  --grid-shape 64 \
  --bounds -5:5 \
  --t-end 5.0 \
  --ic gaussian \
  --no-plot
