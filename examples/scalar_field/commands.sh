#!/usr/bin/env bash
# CLI equivalents for the Klein-Gordon 1+1D example
# See also: klein_gordon.wls (manual derivation), kg_from_lagrangian.py (Python simulation)

set -euo pipefail

# Derive equations from Lagrangian (requires wolframscript)
tg derive theory.toml

# Or preview the generated Wolfram script without running:
# tg derive theory.toml --dry-run

# Inspect the equation system
tg inspect ../data/klein_gordon_1d.json

# Run simulation (Gaussian pulse, matches kg_from_lagrangian.py defaults)
tg simulate ../data/klein_gordon_1d.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-width 5.0 \
  --t-end 30.0 \
  --dt 0.01
