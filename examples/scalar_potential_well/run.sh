#!/usr/bin/env bash
# Scalar Field in External Potential — Background Field Example
#
# Demonstrates the [[background_fields]] feature with two variants:
#   1. Uniform potential: KG equation with constant mass (energy-conserving)
#   2. Localized potential: position-dependent mass via UnitStep
#
# Running:
#   cd examples/scalar_potential_well && uv run bash run.sh

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
  --scheme scipy \
  --output ../data/scalar_potential_well_output

# Step 4: Visualize
tidal plot ../data/scalar_potential_well_output --type heatmap --output ../data/pw_heatmap.png --quiet
tidal plot ../data/scalar_potential_well_output --type profile --output ../data/pw_profile.png --quiet
tidal plot ../data/scalar_potential_well_output --type compare --output ../data/pw_compare.png --quiet
tidal plot ../data/scalar_potential_well_output --type amplitude --output ../data/pw_amplitude.png --quiet

# NOTE: Energy measurement will raise ValueError for the localized variant
# because virial formula doesn't support position-dependent mass terms.
# The simulation itself runs correctly.
echo "Simulation complete (energy not measured — position-dependent mass)"
