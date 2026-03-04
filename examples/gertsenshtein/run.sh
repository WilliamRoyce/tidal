#!/usr/bin/env bash
# Gertsenshtein Effect — Derive → Inspect → Simulate → Measure
#
# Physics: Graviton-photon conversion in uniform background magnetic field B0.
# The Einstein-Maxwell Lagrangian is linearized around flat spacetime + B0,
# producing coupled wave equations for h_+ (graviton) and a_y (photon).
#
# Initial condition: Gaussian wave packet on the graviton (h_+) field,
# observe conversion to photon (a_y) via the Gertsenshtein effect.
#
# Validation: P(graviton -> photon) = sin^2(kappa * B0 * D / (4*sqrt(pi)))
# (Palessandro & Rothman 2023, arXiv:2301.02072, Eq. 26)
#
# Running:
#   cd examples/gertsenshtein && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein Effect: Graviton-Photon Conversion ==="
echo ""

# Step 1: Derive coupled equations from Einstein-Maxwell Lagrangian
echo "--- Step 1: Derive ---"
tidal derive theory.toml
echo ""

# Step 2: Inspect the derived JSON
echo "--- Step 2: Inspect ---"
tidal inspect ../data/gertsenshtein.json
echo ""

# Step 3: Simulate — graviton wave packet propagating through B0 region
# The plane-wave reduction gives 1+1D system along z.
# IC: Gaussian on h_+ (graviton), zero on a_y (photon).
# Expect gradual conversion h_+ -> a_y over the propagation distance.
echo "--- Step 3: Simulate ---"
tidal simulate ../data/gertsenshtein.json \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-amplitude 0.1 \
  --ic-width 3.0 \
  --ic-center 20.0 \
  --t-end 40.0 \
  --param kappa=1.0 B0=0.3 \
  --output ../data/gertsenshtein_output

echo ""

# Step 4: Measure energy conservation and conversion
echo "--- Step 4: Measure ---"
tidal measure ../data/gertsenshtein_output --what conservation --param kappa=1.0 B0=0.3
tidal measure ../data/gertsenshtein_output --what energy --param kappa=1.0 B0=0.3
echo ""

echo "=== Done ==="
