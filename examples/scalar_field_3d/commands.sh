#!/usr/bin/env bash
# CLI equivalents for the Klein-Gordon 3+1D example
# See also: klein_gordon_3d.wls (manual derivation), kg_3d_simulation.py (Python simulation)

set -euo pipefail

# Derive equations from Lagrangian (requires wolframscript)
tg derive theory.toml

# Inspect the equation system
tg inspect ../data/klein_gordon_3d.json

# Run simulation (3D Gaussian pulse at domain center)
tg simulate ../data/klein_gordon_3d.json \
  --param m2=1.0 \
  --grid-shape 32 \
  --bounds 0:20 \
  --periodic \
  --ic gaussian \
  --ic-width 2.0 \
  --t-end 8.0 \
  --dt 0.05
