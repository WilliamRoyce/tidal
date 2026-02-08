#!/usr/bin/env bash
# CLI equivalents for the Polar Klein-Gordon 2+1D example
# See also: polar_kg.wls (manual derivation), polar_kg_simulation.py (Python simulation)
#
# NOTE: This example uses polar coordinates (r, theta) with a coordinate-dependent
# metric. The derive step requires the manual .wls script (no TOML equivalent).
# The simulation uses mixed boundary conditions (Neumann in r, periodic in theta)
# which the CLI does not currently support. The Python script is recommended for
# full-fidelity simulation.

set -euo pipefail

# Derive equations (manual .wls only — curvilinear coordinates)
# tg derive polar_kg.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/polar_kg.json

# Approximate simulation with periodic BCs (loses radial boundary accuracy)
# The Python script uses Neumann BCs in r and periodic in theta
tg simulate ../data/polar_kg.json \
  --param polm2=0.5 \
  --grid-shape 128 \
  --bounds 0.5:10,0:6.283185 \
  --periodic \
  --ic gaussian \
  --ic-center 3.0,3.14159 \
  --ic-width 0.5 \
  --t-end 8.0 \
  --dt 0.005 \
  --scheme scipy
