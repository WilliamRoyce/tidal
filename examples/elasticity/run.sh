#!/usr/bin/env bash
# CLI equivalents for the Elasticity (Navier-Cauchy) 2+1D example
# See also: elasticity_from_lagrangian.py (Python simulation)
#
# The elasticity Lagrangian uses component-derivative notation (CD[{idx, -chart}])
# because the Lame parameters produce anisotropic spatial coefficients.
#
# To run manually:  cd examples/elasticity

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from TOML config (uses component-derivative -chart placeholder)
tidal derive theory.toml

# Inspect the equation system
tidal inspect ../data/navier_cauchy_2d.json

# Run simulation (Gaussian pulse in ux displacement)
# Parameters rho, lambda, mu are baked into the JSON as numeric coefficients
tidal simulate ../data/navier_cauchy_2d.json \
  --grid-shape 64 \
  --bounds 0:10 \
  --periodic \
  --ic gaussian \
  --ic-component ux_0 \
  --ic-width 1.0 \
  --t-end 3.0 \
  --dt 0.005
