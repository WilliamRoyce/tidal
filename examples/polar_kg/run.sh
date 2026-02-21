#!/usr/bin/env bash
# Polar Klein-Gordon 2+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: Klein-Gordon in polar coordinates (r, θ). The coordinate-dependent
# metric (ds² = -dt² + dr² + r²dθ²) produces position-dependent coefficients
# and Christoffel corrections. Neumann BC in r, periodic in θ.
#
# Running this script:
#   cd examples/polar_kg && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/polar_kg.json
#   tidal simulate ../data/polar_kg.json --param polm2=0.5 \
#     --grid-shape 128 --bounds 0.5:10,0:6.283185 --bc neumann,periodic \
#     --ic formula --ic-formula "np.exp(-(x - 3.0)**2 / 0.5)" \
#     --t-end 8.0 --scheme scipy --output ../data/polar_kg_output
#   tidal plot ../data/polar_kg_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/polar_kg.json

# Run simulation (Gaussian ring at r=3, Neumann in r, periodic in theta)
# Coordinates: x=r, y=theta
tidal simulate ../data/polar_kg.json \
  --param polm2=0.5 \
  --grid-shape 128 \
  --bounds 0.5:10,0:6.283185 \
  --bc neumann,periodic \
  --ic formula \
  --ic-formula "np.exp(-(x - 3.0)**2 / 0.5)" \
  --t-end 8.0 \
  --scheme scipy \
  --output ../data/polar_kg_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/polar_kg_output --type snapshot --time-index 0 --quiet
tidal plot ../data/polar_kg_output --type snapshot --time-index -1 --quiet
tidal plot ../data/polar_kg_output --type amplitude --quiet
