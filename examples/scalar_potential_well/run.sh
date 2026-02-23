#!/usr/bin/env bash
# Scalar Field in External Potential 1+1D — Background Field Example
#
# Physics: Demonstrates the [[background_fields]] feature with a localized
# potential well via UnitStep. A Gaussian wave packet scatters off the
# potential boundary, producing reflected and transmitted components.
# Position-dependent mass from the external potential V(x).
#
# Running this script:
#   cd examples/scalar_potential_well && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/scalar_potential_well.json
#   tidal simulate ../data/scalar_potential_well.json --param V0=4.0 \
#     --grid-shape 256 --bounds 0:100 --periodic --ic gaussian \
#     --ic-component phi_0 --ic-center 15.0 --ic-width 5.0 \
#     --t-end 20.0 --output ../data/scalar_potential_well_output
#   tidal plot ../data/scalar_potential_well_output --type heatmap --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive
tidal derive theory.toml

# Step 2: Inspect
tidal inspect ../data/scalar_potential_well.json

# Step 3: Simulate — Gaussian pulse scatters off potential boundary
tidal simulate ../data/scalar_potential_well.json \
  --param V0=4.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-component phi_0 \
  --ic-center 15.0 \
  --ic-width 5.0 \
  --t-end 20.0 \
  --output ../data/scalar_potential_well_output

# Step 4: Visualize (plots saved into the simulation output directory)
tidal plot ../data/scalar_potential_well_output --type heatmap --quiet
tidal plot ../data/scalar_potential_well_output --type profile --quiet
tidal plot ../data/scalar_potential_well_output --type compare --quiet
tidal plot ../data/scalar_potential_well_output --type amplitude --quiet

# NOTE: Energy measurement will raise ValueError for the localized variant
# because virial formula doesn't support position-dependent mass terms.
# The simulation itself runs correctly.
echo "Simulation complete (energy not measured — position-dependent mass)"
