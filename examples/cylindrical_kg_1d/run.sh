#!/usr/bin/env bash
# Cylindrical Klein-Gordon (axial plane wave) — Reduced 1+1D pipeline
#
# Physics: Klein-Gordon in cylindrical coordinates, plane-wave reduction
# along z-axis (axial propagation). Since the z-axis in cylindrical
# coordinates is flat (g_zz = 1, no Christoffels), this reduces to
# a simple 1D massive wave equation: d2t(phi) = d2z(phi) - m2*phi.
# No volume element (flat axis), no position-dependent coefficients.
#
# This tests that reduction correctly eliminates all curved-coordinate
# complexity when the propagation axis is flat.
#
# Running this script:
#   cd examples/cylindrical_kg_1d && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/cylindrical_kg_1d.json
#   tidal simulate ../data/cylindrical_kg_1d.json \
#     --grid-shape 128 --bounds -5:5 --bc neumann \
#     --ic gaussian --ic-width 0.5 --t-end 4.0
#   tidal measure output/ --what conservation

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system (should show 1+1D, 1 field)
tidal inspect ../data/cylindrical_kg_1d.json

# Run 1D simulation (Gaussian pulse on flat z-axis)
tidal simulate ../data/cylindrical_kg_1d.json \
  --grid-shape 128 \
  --bounds -5:5 \
  --bc neumann \
  --ic gaussian \
  --ic-width 0.5 \
  --t-end 4.0

# Check energy conservation
tidal measure output/ --what conservation
