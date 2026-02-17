#!/usr/bin/env bash
# CLI equivalents for the Klein-Gordon 1+1D example
# See also: kg_from_lagrangian.py (Python simulation)
#
# To run manually:  cd examples/scalar_field && tidal derive theory.toml

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Or preview the generated Wolfram script without running:
# tidal derive theory.toml --dry-run

# Inspect the equation system
tidal inspect ../data/klein_gordon_1d.json

# Run simulation (Gaussian pulse, matches kg_from_lagrangian.py defaults)
tidal simulate ../data/klein_gordon_1d.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-width 5.0 \
  --t-end 30.0 \
  --dt 0.01 \
  --snapshots 0.1
