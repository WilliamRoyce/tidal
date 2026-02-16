#!/usr/bin/env bash
# CLI equivalents for the Electrostatics 2D example
#
# Derives the free scalar field wave equation from a Lagrangian,
# then demonstrates constraint-mode (Laplace equation) solving.
#
# To run manually:  cd examples/electrostatics

set -euo pipefail
cd "$(dirname "$0")"

# Derive wave equation from TOML config (Lagrangian-derived)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/electrostatics_laplace.json

# Simulate: wave equation with Gaussian initial condition
tidal simulate ../data/electrostatics_laplace.json \
  --grid-shape 64 \
  --bounds -5:5 \
  --t-end 5.0 \
  --ic gaussian \
  --no-plot
