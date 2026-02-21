#!/usr/bin/env bash
# Elasticity (Navier-Cauchy) 2+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: The elasticity Lagrangian uses component-derivative notation
# (CD[{idx, -chart}]) because the Lamé parameters produce anisotropic spatial
# coefficients. Two displacement fields (ux, uy) are coupled via shear.
#
# Running this script:
#   cd examples/elasticity && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/navier_cauchy_2d.json
#   tidal simulate ../data/navier_cauchy_2d.json --grid-shape 64 --bounds 0:10 \
#     --periodic --ic gaussian --ic-component ux_0 --ic-width 1.0 --t-end 3.0 \
#     --scheme scipy --output ../data/elasticity_output
#   tidal plot ../data/elasticity_output --type amplitude --quiet

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
  --scheme scipy \
  --output ../data/elasticity_output

# Visualize results — snapshots of each displacement component
tidal plot ../data/elasticity_output --type snapshot --field ux_0 --time-index -1 --quiet
tidal plot ../data/elasticity_output --type snapshot --field uy_0 --time-index -1 --quiet
tidal plot ../data/elasticity_output --type amplitude --quiet
