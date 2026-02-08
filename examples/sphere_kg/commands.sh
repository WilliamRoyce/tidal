#!/usr/bin/env bash
# CLI equivalents for the Sphere Klein-Gordon 2+1D example (stereographic projection)
# See also: sphere_kg.wls (manual derivation), sphere_kg_simulation.py (Python simulation)
#
# NOTE: The derive step uses stereographic projection coordinates with
# position-dependent metric. The manual .wls script is required for derivation.
# The simulate step works fully via CLI (all periodic BCs, 2D Cartesian grid).
#
# To run manually:  cd examples/sphere_kg

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations (manual .wls only — stereographic projection metric)
# tg derive sphere_kg.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/sphere_kg.json

# Run simulation (2D Gaussian on stereographic plane)
# Position-dependent wave speed from conformal factor Omega(r) = 2R^2/(R^2 + r^2)
tg simulate ../data/sphere_kg.json \
  --param sR=2.0 \
  --grid-shape 128 \
  --bounds -8:8 \
  --periodic \
  --ic gaussian \
  --ic-width 0.8 \
  --t-end 10.0 \
  --dt 0.005
