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
#   tidal simulate ../data/navier_cauchy_2d.json --grid-shape 128 --bounds 0:10 \
#     --periodic --ic gaussian --ic-component ux_0 --ic-width 1.0 --t-end 3.0 \
#     --output ../data/elasticity_output
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
  --grid-shape 128 \
  --bounds 0:10 \
  --periodic \
  --ic gaussian \
  --ic-component ux_0 \
  --ic-width 1.0 \
  --t-end 3.0 \
  --output ../data/elasticity_output

# Visualize results — snapshots of each displacement component
tidal plot ../data/elasticity_output --type snapshot --field ux_0 --time-index -1 --quiet
tidal plot ../data/elasticity_output --type snapshot --field uy_0 --time-index -1 --quiet
tidal plot ../data/elasticity_output --type amplitude --quiet
tidal plot ../data/elasticity_output --type hamiltonian --quiet
tidal plot ../data/elasticity_output --type conservation --quiet

echo ""
echo "=== Run 1 (compression wave in ux) complete ==="

# ============================================================================
# Run 2: Shear wave excitation (using --ic-field for velocity)
# ============================================================================
# Run 1 excites longitudinal compression (displacement in ux). Elasticity's
# interesting physics is the coupling between longitudinal and shear modes.
# Using --ic-field with a velocity formula on uy_0, we launch a transverse
# (shear) wave — zero displacement, nonzero velocity — producing a
# qualitatively different wave pattern.

tidal simulate ../data/navier_cauchy_2d.json \
  --grid-shape 128 \
  --bounds 0:10 \
  --periodic \
  --scheme leapfrog \
  --ic zero \
  --ic-field 'uy_0:velocity:exp(-(x - 5)**2 / 0.5) * cos(4 * y)' \
  --t-end 3.0 \
  --output ../data/elasticity_shear

# Compare compression (Run 1) vs shear (Run 2) patterns
tidal plot ../data/elasticity_shear --type snapshot --field ux_0 --time-index -1 \
  --title "ux (shear excitation)" --quiet
tidal plot ../data/elasticity_shear --type snapshot --field uy_0 --time-index -1 \
  --title "uy (shear excitation)" --quiet
tidal plot ../data/elasticity_shear --type amplitude \
  --title "Shear wave amplitudes" --quiet
tidal plot ../data/elasticity_shear --type conservation --quiet

echo ""
echo "=== Run 2 (shear wave) complete ==="
echo ""
echo "Compare compression vs shear wave patterns:"
echo "  Compression (ux pulse):  ../data/elasticity_output/"
echo "  Shear (uy velocity):     ../data/elasticity_shear/"
