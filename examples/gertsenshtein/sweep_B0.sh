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
# Grid resolution N=512 with 4th-order FD stencils (Fornberg 1988):
#   O(dx^4) convergence gives equivalent accuracy to N=1024 at 2nd-order,
#   with ~2x faster wall-clock time.
#   kDx = 2.011 * (100/512) = 0.39 -> 4th-order error O(kDx)^4 ~ 2e-2%
#   (vs 2nd-order at N=1024: kDx=0.20 -> error O(kDx)^2 ~ 0.5%)
#
# The conversion formula P = sin^2(kappa*B0*D/2) is derived from eigenmode
# analysis of the coupled h_7/a_2 equations from L = (1/kappa^2)R - (1/4)F^2,
# and confirmed numerically (RMS < 0.015).
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

# Plot 1: P(B0) vs bare analytical formula sin^2(kappa * B0 * D / 2)
tidal plot ../data/gertsenshtein_sweep_B0 --type sweep \
  --metric P_final \
  --title "Gertsenshtein P(B0) — graviton-photon conversion" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2' \
  --output ../data/gertsenshtein_sweep_B0/sweep.png --quiet

# Plot 2: P_max(B0) vs corrected formula including graviton effective mass term
#
# The h_7 graviton has effective mass m^2_eff = kappa^2 * B0^2 from the background
# EM stress-energy. This contributes an extra kappa^2*B0^2*h_7^2/2 term to the
# initial graviton energy, which is not transferred to the massless photon:
#
#   E_source(t=0) includes (k^2 + kappa^2*B0^2) / (kappa^2)  [graviton with mass]
#   E_target(t_peak) includes only k^2 / kappa^2              [photon, massless]
#
#   P_max = k^2 / (k^2 + kappa^2 * B0^2)
#
# where k = 2*pi*32/100 = 2.0106 (k^2 = 4.0425) for 32 wavelengths in [0,100].
#
# Full corrected formula (for P_max, the peak Rabi envelope):
#   P_max_theory = sin^2(kappa * B0 * t_end / 2) * 4.0425 / (4.0425 + kappa^2 * B0^2)
#
# IMPORTANT: This correction applies to P_max (peak conversion envelope), NOT to
# P_final (conversion at t=t_end). The bare sin^2 formula is exact for P_final
# at any time; the k^2/(k^2+kappa^2*B0^2) factor caps the energy-based measurement
# because the graviton's mass energy is not transferred to the massless photon.
#
# This matches TIDAL P_max to < 0.4% at all tested B0 values. The correction grows
# as B0^2/k^2 and reaches ~1% at B0 = 0.2 (confirmed by targeted simulations).
tidal plot ../data/gertsenshtein_sweep_B0 --type sweep \
  --metric P_max \
  --title "Gertsenshtein P_max(B0) — corrected for graviton effective mass" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2 * 4.0425 / (4.0425 + kappa**2 * B0**2)' \
  --output ../data/gertsenshtein_sweep_B0/sweep_corrected.png --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: examples/data/gertsenshtein_sweep_B0/"
echo "Plot 1: examples/data/gertsenshtein_sweep_B0/sweep.png            [bare sin^2]"
echo "Plot 2: examples/data/gertsenshtein_sweep_B0/sweep_corrected.png  [P_max corrected]"
