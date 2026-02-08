#!/usr/bin/env bash
# CLI equivalents for the Cylindrical Klein-Gordon 3+1D example
# See also: cylindrical_kg.wls (manual derivation), cylindrical_kg_simulation.py (Python simulation)
#
# NOTE: This example uses cylindrical coordinates (r, theta, z) with a
# coordinate-dependent metric. The derive step requires the manual .wls script.
# The simulation uses mixed BCs (Neumann in r and z, periodic in theta)
# which the CLI does not support. The Python script is recommended.

set -euo pipefail

# Derive equations (manual .wls only — curvilinear coordinates)
# tg derive cylindrical_kg.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/cylindrical_kg.json

# Approximate simulation with periodic BCs (loses radial/axial boundary accuracy)
tg simulate ../data/cylindrical_kg.json \
  --grid-shape 48 \
  --bounds 0.5:8,0:6.283185,-5:5 \
  --periodic \
  --ic gaussian \
  --ic-center 3.0,3.14159,0.0 \
  --ic-width 0.6 \
  --t-end 4.0 \
  --dt 0.01
