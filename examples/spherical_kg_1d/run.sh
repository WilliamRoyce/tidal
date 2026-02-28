#!/usr/bin/env bash
# Spherical Klein-Gordon (radial plane wave) — Reduced 1+1D pipeline
#
# Physics: Klein-Gordon in spherical coordinates, plane-wave reduction
# along x(=r) (radial propagation, d/d(theta) = d/d(phi) = 0).
# The reduced EOM: d2t(phi) = d2r(phi) + (2/r)*dr(phi) - m2*phi.
# Volume element: r^2 (from factored sqrt|det(g_spatial)|).
#
# Energy conservation: E = integral[ (v^2 + (dr phi)^2 + m2*phi^2) r^2 dr ]
# is exactly conserved (verified analytically: boundary terms vanish).
#
# This tests curved-coordinate reduction with:
#   - Position-dependent coefficient (2/r from Christoffel)
#   - Non-trivial volume element (r^2)
#   - All angular operators eliminated
#
# Running this script:
#   cd examples/spherical_kg_1d && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/spherical_kg_1d.json
#   tidal simulate ../data/spherical_kg_1d.json \
#     --grid-shape 128 --bounds 0.5:10 --bc neumann \
#     --ic gaussian --ic-width 0.5 --ic-center 4.0 --t-end 5.0
#   tidal measure output/ --what conservation

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system (should show 1+1D with volume_element)
tidal inspect ../data/spherical_kg_1d.json

# Run 1D simulation (Gaussian shell at r=4, Neumann BCs)
# Bounds start at 0.5 to avoid r=0 singularity in 2/r coefficient
tidal simulate ../data/spherical_kg_1d.json \
  --grid-shape 128 \
  --bounds 0.5:10 \
  --bc neumann \
  --ic gaussian \
  --ic-width 0.5 \
  --ic-center 4.0 \
  --t-end 5.0

# Check energy conservation (should be < 1% for 128-point grid)
tidal measure output/ --what conservation
