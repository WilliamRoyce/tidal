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
#     --t-end 8.0 --output ../data/polar_kg_output
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
  --output ../data/polar_kg_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/polar_kg_output --type snapshot --time-index 0 --quiet
tidal plot ../data/polar_kg_output --type snapshot --time-index -1 --quiet
tidal plot ../data/polar_kg_output --type amplitude --quiet

echo ""
echo "=== Run 1 (static ring) complete ==="

# ============================================================================
# Run 2: Expanding ring (outward radial velocity)
# ============================================================================
# Run 1 creates a static Gaussian ring at r=3 which collapses inward and
# expands outward simultaneously. Using --ic-formula-velocity with a positive
# radial velocity, we create a purely expanding ring — a more physically
# natural IC on curved manifolds.
#
# In polar coords, x=r and positive velocity means outward expansion.

tidal simulate ../data/polar_kg.json \
  --param polm2=0.5 \
  --grid-shape 128 \
  --bounds 0.5:10,0:6.283185 \
  --bc neumann,periodic \
  --ic formula \
  --ic-formula 'exp(-(x - 3.0)**2 / 0.5)' \
  --ic-formula-velocity '0.5 * exp(-(x - 3.0)**2 / 0.5)' \
  --t-end 8.0 \
  --output ../data/polar_kg_expanding

# Compare: static (symmetric collapse+expansion) vs expanding (outward only)
tidal plot ../data/polar_kg_expanding --type snapshot --time-index 0 \
  --title "Expanding ring t=0" --quiet
tidal plot ../data/polar_kg_expanding --type snapshot --time-index -1 \
  --title "Expanding ring t=final" --quiet
tidal plot ../data/polar_kg_expanding --type amplitude \
  --title "Expanding ring amplitude" --quiet

echo ""
echo "=== Run 2 (expanding ring) complete ==="
echo ""
echo "Compare static vs expanding ring:"
echo "  Static (no velocity):    ../data/polar_kg_output/"
echo "  Expanding (v = +0.5):    ../data/polar_kg_expanding/"
