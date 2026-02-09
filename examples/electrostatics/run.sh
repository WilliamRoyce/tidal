#!/usr/bin/env bash
# CLI equivalents for the Electrostatics 2D example
# See also: electrostatics.wls (manual derivation), electrostatics_simulation.py (Python simulation)
#
# Electrostatics solves the Poisson equation (no time evolution).
# Use --mode constraint for single-step constraint solving.
#
# To run manually:  cd examples/electrostatics

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations (Poisson equation, .wls pass-through to wolframscript)
tidal derive electrostatics.wls

# Inspect the equation system
tidal inspect ../data/electrostatics_2d.json

# Solve constraint (Poisson equation for phi given charge density rho)
# Gaussian charge density at origin, Dirichlet BCs (phi=0 at boundary)
tidal simulate ../data/electrostatics_2d.json \
  --mode constraint \
  --grid-shape 64 \
  --bounds -5:5 \
  --bc dirichlet \
  --ic formula \
  --ic-formula "-np.exp(-(x**2 + y**2) / 0.5)" \
  --ic-component rho
