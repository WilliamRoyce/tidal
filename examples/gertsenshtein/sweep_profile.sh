#!/usr/bin/env bash
# Gertsenshtein Localized B-field Sweep — Validation of Boccaletti Formula
#
# Sweeps Bpeak (peak field strength) and R (Gaussian half-width) and measures
# the graviton-photon conversion probability P(graviton -> photon) at each point.
#
# NOTE: constant "Bpeak" (not "B0_peak") — underscores are Mathematica pattern syntax.
#
# Analytical target (Boccaletti 1970, weak-field Gaussian):
#   P = sin^2(kappa * Bpeak * R * sqrt(pi/2))
#     = sin^2(kappa * integral(B_x(z), -inf, +inf) / 2)
#
# where B_x(z) = Bpeak * exp(-z^2 / (2*R^2)) and integral = Bpeak * R * sqrt(2*pi).
#
# Parameter ranges:
#   Bpeak in [0.02, 0.20] (8 points) — scan from weak to moderate coupling
#   R in [2.0, 15.0] (6 points)      — scan from narrow to wide profiles
#
# For kappa=1, the argument kappa*Bpeak*R*sqrt(pi/2) ranges from
#   min: 1.0 * 0.02 * 2.0 * 1.253 = 0.050  (P ≈ 0.0025 — deep weak-field)
#   max: 1.0 * 0.20 * 15.0 * 1.253 = 3.758  (P oscillates — beyond first maximum)
# This covers the full weak-field → first Rabi maximum regime.
#
# Simulation setup:
#   - Non-periodic domain [-100, 100] (Neumann BCs) — periodic is incorrect for Abar_y
#   - Gaussian wavepacket at z = -50, σ=5, k = 2.0, t_end = 120
#   - t_end=120 ensures exactly one passage through B-field (wavepacket hits right wall at t≈150)
#   - Conversion measured at peak P(t) during single pass through B-field
#
# Note: Use --what conversion (not --what energy) — energy measurement is unreliable
# with position-dependent coefficients (known limitation of current implementation).
#
# Ref: Boccaletti et al. (1970, Nuovo Cimento 70B:129)
#      Domcke, Garcia-Cely & Lee (2025, arXiv:2507.16609)
#
# Running:
#   cd examples/gertsenshtein && bash sweep_profile.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein Localized B-field Profile Sweep ==="
echo ""

tidal sweep ../data/gertsenshtein_localized.json \
  --sweep "Bpeak=0.02:0.20:8" \
  --sweep "R=2.0:15.0:6" \
  --measure conversion \
  --source h_7 --target a_2 \
  --grid-shape 1024 \
  "--bounds=-100:100" \
  --ic gaussian \
  --ic-wavevector 2.0 \
  --ic-amplitude 0.1 \
  --ic-width 5.0 \
  "--ic-center=-50.0" \
  --ic-component h_7 \
  --t-end 120.0 \
  --param kappa=1.0 \
  --output ../data/gertsenshtein_sweep_profile

echo ""

# Plot sweep results: P(B0_peak, R) vs Boccaletti formula
tidal plot ../data/gertsenshtein_sweep_profile --type sweep \
  --metric P_final \
  --title "Gertsenshtein P(Bpeak, R) — Boccaletti formula validation" \
  --output ../data/gertsenshtein_sweep_profile/sweep.png --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: examples/data/gertsenshtein_sweep_profile/"
echo "Plot:    examples/data/gertsenshtein_sweep_profile/sweep.png"
