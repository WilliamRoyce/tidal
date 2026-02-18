#!/usr/bin/env bash
# CLI equivalents for the Chern-Simons 2+1D example
# See also: chern_simons_simulation.py (Python simulation)
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
  --scheme scipy
