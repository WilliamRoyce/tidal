#!/usr/bin/env bash
# CLI equivalents for the Massive 3-Form 3+1D example
# See also: massive_3form.wls (manual derivation), simulate_massive_3form.py (Python simulation)

set -euo pipefail

# Derive equations from Lagrangian (requires wolframscript)
tg derive theory.toml

# Inspect the equation system
tg inspect ../data/massive_3form.json

# Run simulation (Gaussian pulse in C_0, other components start at zero)
# Antisymmetric rank-3 tensor: 64 components reduce to 4 independent (C_0..C_3)
tg simulate ../data/massive_3form.json \
  --param m2=1.0 \
  --grid-shape 16 \
  --bounds 0:10 \
  --periodic \
  --ic gaussian \
  --ic-component C_0 \
  --ic-width 1.5 \
  --t-end 5.0 \
  --dt 0.05
