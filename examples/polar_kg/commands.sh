#!/usr/bin/env bash
# CLI equivalents for the Polar Klein-Gordon 2+1D example
# See also: polar_kg.wls (manual derivation), polar_kg_simulation.py (Python simulation)
#
# NOTE: The derive step uses polar coordinates (r, theta) with a coordinate-dependent
# metric. The manual .wls script is required for derivation (no TOML equivalent).
# The simulation works via CLI with --bc and --ic formula flags.
#
# To run manually:  cd examples/polar_kg

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations (manual .wls only — curvilinear coordinates)
# tg derive polar_kg.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/polar_kg.json

# Run simulation (Gaussian ring at r=3, Neumann in r, periodic in theta)
# Coordinates: x=r, y=theta
tg simulate ../data/polar_kg.json \
  --param polm2=0.5 \
  --grid-shape 128 \
  --bounds 0.5:10,0:6.283185 \
  --bc neumann,periodic \
  --ic formula \
  --ic-formula "np.exp(-(x - 3.0)**2 / 0.5)" \
  --t-end 8.0 \
  --dt 0.005 \
  --scheme scipy
