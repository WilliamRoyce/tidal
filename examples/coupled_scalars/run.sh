#!/usr/bin/env bash
# CLI equivalents for the Coupled Scalars 1+1D example
# See also: coupled_scalars.wls (manual derivation), coupled_from_lagrangian.py (Python simulation)
#
# To run manually:  cd examples/coupled_scalars && tg derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tg derive theory.toml

# Inspect the equation system
tg inspect ../data/coupled_scalars.json

# Run simulation (off-center Gaussian in phi only, chi starts at zero)
# Coupling transfers energy from phi to chi over time
tg simulate ../data/coupled_scalars.json \
  --param mPhi2=1.0 --param mChi2=4.0 --param gCpl=0.5 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component phi_0 \
  --ic-center 30.0 \
  --ic-width 5.0 \
  --t-end 20.0 \
  --dt 0.01
