#!/usr/bin/env bash
# Gertsenshtein B0 Sweep — Validation of P = sin^2(kappa * B0 * D / 2)
#
# Sweeps B0 (background magnetic field strength) and measures the conversion
# probability P(graviton -> photon) at each value. The results should trace
# the analytical Rabi oscillation curve.
#
# IC: monochromatic plane wave cos(k*z) on h_7 (graviton h_+)
# k = 2*pi*32/100 ~ 2.011 (32 wavelengths in domain [0,100])
#
# This ensures k >> B0*kappa (massless dispersion limit) across the sweep
# range B0 in [0.005, 0.25]:
#   k/(B0_max * kappa) = 2.011/0.25 ~ 8 >> 1
#   Mass correction O(B0^2/k^2) < 1.6% at B0_max
#
# With D = c * t_end = 50, kappa = 1:
#   P(B0) = sin^2(B0 * 50 / 2) = sin^2(25 * B0)
#   First maximum at B0 = pi/50 ~ 0.063
#   Period: delta(B0) = pi/25 ~ 0.126
#   Sweep covers ~2 full oscillation cycles
#
# The conversion formula P = sin^2(kappa*B0*D/2) is derived from eigenmode
# analysis of the coupled h_7/a_2 equations from L = (1/kappa^2)R - (1/4)F^2,
# and confirmed numerically (RMS < 0.02).
#
# Ref: Gertsenshtein (1962), JETP 14:84
#
# Running:
#   cd examples/gertsenshtein && bash sweep_B0.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein B0 Coupling Sweep ==="
echo ""

tidal sweep ../data/gertsenshtein.json \
  --sweep "B0=0.005:0.25:40" \
  --measure conversion \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-wavevector 2.0106 \
  --ic-amplitude 0.1 \
  --ic-component h_7 \
  --t-end 50.0 \
  --param kappa=1.0 \
  --output ../data/gertsenshtein_sweep_B0

echo ""

# Plot sweep results: P(B0) vs analytical
tidal plot ../data/gertsenshtein_sweep_B0 --type sweep \
  --metric P_final \
  --title "Gertsenshtein P(B0) — graviton-photon conversion" \
  --output ../data/gertsenshtein_sweep_B0/sweep.png --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: examples/data/gertsenshtein_sweep_B0/"
echo "Plot:   examples/data/gertsenshtein_sweep_B0/sweep.png"
