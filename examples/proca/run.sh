#!/usr/bin/env bash
# CLI equivalents for the Proca (massive vector) 1+1D example
# See also: proca_simulation.py (Python simulation)
#
# To run manually:  cd examples/proca && tidal derive theory.toml

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
  --dt 0.005
