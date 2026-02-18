#!/usr/bin/env bash
# CLI equivalents for the Cylindrical Klein-Gordon 3+1D example
# See also: cylindrical_kg_simulation.py (Python simulation)
#
# NOTE: The derive step uses cylindrical coordinates (r, theta, z) with a
# coordinate-dependent metric.
# The simulation works via CLI with --bc and --ic formula flags.
#
# To run manually:  cd examples/cylindrical_kg

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/cylindrical_kg.json

# Run simulation (Gaussian ring in r-z plane, Neumann in r and z, periodic in theta)
# Coordinates: x=r, y=theta, z=z
tidal simulate ../data/cylindrical_kg.json \
  --param cylm2=0.5 \
  --grid-shape 48 \
  --bounds 0.5:8,0:6.283185,-5:5 \
  --bc neumann,periodic,neumann \
  --ic formula \
  --ic-formula "np.exp(-((x - 3.0)**2 / 0.72) - (z**2 / 1.28))" \
  --t-end 4.0 \
  --scheme scipy
