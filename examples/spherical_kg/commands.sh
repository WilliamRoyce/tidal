#!/usr/bin/env bash
# CLI equivalents for the Spherical Klein-Gordon 3+1D example
# See also: spherical_kg.wls (manual derivation), spherical_kg_simulation.py (Python simulation)
#
# NOTE: This example uses spherical coordinates (r, theta, phi) with a
# coordinate-dependent metric and trigonometric coefficient functions.
# The derive step requires the manual .wls script (no TOML equivalent).
# The simulation uses mixed BCs (Neumann in r and theta, periodic in phi)
# which the CLI does not support. The Python script is recommended.

set -euo pipefail

# Derive equations (manual .wls only — curvilinear coordinates)
# tg derive spherical_kg.wls    # pass-through to wolframscript

# Inspect the equation system
tg inspect ../data/spherical_kg.json

# Approximate simulation with periodic BCs (loses boundary accuracy at poles/origin)
tg simulate ../data/spherical_kg.json \
  --param sphm2=0.5 \
  --grid-shape 64 \
  --bounds 0.5:8,0.05:3.09,0:6.283185 \
  --periodic \
  --ic gaussian \
  --ic-center 3.0,1.5708,3.14159 \
  --ic-width 0.6 \
  --t-end 5.0 \
  --dt 0.01
