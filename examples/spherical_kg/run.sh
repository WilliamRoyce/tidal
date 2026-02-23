#!/usr/bin/env bash
# Spherical Klein-Gordon 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: Klein-Gordon in spherical coordinates (r, θ, φ). The metric
# ds² = -dt² + dr² + r²dθ² + r²sin²θ dφ² produces trigonometric
# coefficient functions and Christoffel corrections. Neumann in r and θ,
# periodic in φ.
#
# NOTE: 3D data — tidal plot is for 1D/2D only.
#
# Running this script:
#   cd examples/spherical_kg && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/spherical_kg.json
#   tidal simulate ../data/spherical_kg.json --param spm2=0.5 \
#     --grid-shape 64 --bounds 0.5:8,0.05:3.09,0:6.283185 \
#     --bc neumann,neumann,periodic \
#     --ic formula --ic-formula "np.exp(-(x - 3.0)**2 / 0.72)" \
#     --t-end 5.0

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/spherical_kg.json

# Run simulation (Gaussian shell at r=3, Neumann in r and theta, periodic in phi)
# Coordinates: x=r, y=theta, z=phi
tidal simulate ../data/spherical_kg.json \
  --param spm2=0.5 \
  --grid-shape 64 \
  --bounds 0.5:8,0.05:3.09,0:6.283185 \
  --bc neumann,neumann,periodic \
  --ic formula \
  --ic-formula "np.exp(-(x - 3.0)**2 / 0.72)" \
  --t-end 5.0
