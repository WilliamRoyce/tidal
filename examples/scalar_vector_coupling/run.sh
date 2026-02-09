#!/usr/bin/env bash
# CLI equivalents for the Scalar-Vector Coupling stress test example
# See also: scalar_vector_coupling.wls (manual derivation)
#
# This example stress-tests the pipeline with:
#   - Mixed-rank cross-field coupling (scalar phi + vector A)
#   - Cross-field first_derivative_t and gradient operators
#   - Epsilon tensor + cross-field coupling in same system
#   - 4 symbolic constants (phim2, Am2, kCS, gSV)
#   - 4x4 mass/coupling matrices

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/scalar_vector_coupling.json

# Run simulation (Gaussian IC for phi, periodic BCs)
tidal simulate ../data/scalar_vector_coupling.json \
  --param phim2=1.0 --param Am2=0.5 --param kCS=0.3 --param gSV=0.2 \
  --grid-shape 48 \
  --bounds 0:10,0:10 \
  --bc periodic,periodic \
  --ic gaussian \
  --t-end 5.0 \
  --dt 0.005
