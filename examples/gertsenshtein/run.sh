#!/usr/bin/env bash
# Gertsenshtein Effect — Derive → Inspect → Simulate → Measure
#
# Physics: Graviton-photon conversion in uniform background magnetic field B0.
# The Einstein-Maxwell Lagrangian is linearized around flat spacetime + B0,
# producing coupled wave equations for h_+ (graviton) and a_y (photon).
#
# Gauge choices (following Palessandro & Rothman 2023, arXiv:2301.02072):
#   h: TT gauge (transverse-traceless) — standard for linearized gravity
#   a: Lorenz gauge (d_mu a^mu = 0) — decouples EM into wave equations
# Together these reduce the full tensor+vector system to just h_+ and a_y.
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
# After constraint elimination, h_7 = -h_xx is the plus polarization.
# IC: Gaussian on h_7 (graviton h_+), zero on a_2 (photon A_y).
# Expect gradual conversion h_7 -> a_2 over the propagation distance.
echo "--- Step 3: Simulate ---"
tidal simulate ../data/gertsenshtein.json \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-amplitude 0.1 \
  --ic-width 3.0 \
  --ic-center 20.0 \
  --ic-component h_7 \
  --t-end 40.0 \
  --param kappa=1.0 --param B0=0.3 \
  --output ../data/gertsenshtein_output

echo ""

# Step 4: Measure energy conservation and conversion
echo "--- Step 4: Measure ---"
tidal measure ../data/gertsenshtein_output --what conservation --param kappa=1.0 --param B0=0.3
tidal measure ../data/gertsenshtein_output --what energy --param kappa=1.0 --param B0=0.3
echo ""

echo "=== Done ==="
