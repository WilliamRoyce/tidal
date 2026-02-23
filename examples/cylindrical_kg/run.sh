#!/usr/bin/env bash
# Cylindrical Klein-Gordon 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: Klein-Gordon in cylindrical coordinates (r, θ, z). The metric
# ds² = -dt² + dr² + r²dθ² + dz² produces Christoffel corrections and
# position-dependent coefficients. Neumann in r and z, periodic in θ.
#
# NOTE: 3D data — tidal plot is for 1D/2D only.
#
# Running this script:
#   cd examples/cylindrical_kg && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/cylindrical_kg.json
#   tidal simulate ../data/cylindrical_kg.json --param cylm2=0.5 \
#     --grid-shape 48 --bounds 0.5:8,0:6.283185,-5:5 \
#     --bc neumann,periodic,neumann \
#     --ic formula --ic-formula "np.exp(-((x - 3.0)**2 / 0.72) - (z**2 / 1.28))" \
#     --t-end 4.0

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
  --t-end 4.0
