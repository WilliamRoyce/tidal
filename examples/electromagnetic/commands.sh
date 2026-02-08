#!/usr/bin/env bash
# CLI equivalents for the Electromagnetic 1+1D example
# See also: em_lagrangian_1d.wls (manual derivation), em_from_lagrangian.py (Python simulation)
#
# To run manually:  cd examples/electromagnetic && tg derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tg derive theory.toml

# Inspect the equation system
tg inspect ../data/em_1d.json

# Run simulation (Gaussian pulse in A_1 component)
# A_0 is a constraint (time_order=0), only A_1 is dynamical
tg simulate ../data/em_1d.json \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component A_1 \
  --ic-width 5.0 \
  --t-end 25.0 \
  --dt 0.01
