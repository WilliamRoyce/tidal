#!/usr/bin/env bash
# Gertsenshtein Effect (Localized B-field) — Derive → Inspect → Simulate → Measure
#
# Physics: Graviton-photon conversion in a Gaussian background magnetic field.
# B_x(z) = Bpeak * exp(-z^2 / (2*R^2)), localized around z=0.
# Gauge potential: Abar_y = -Bpeak * R * sqrt(pi/2) * erf(z / (sqrt(2)*R))
#
# NOTE: constant "Bpeak" (not "B0_peak") — underscores are Mathematica pattern syntax.
#
# An incident graviton wavepacket (h_7) starts at z=-50, propagates through
# the magnetized region at z=0, and partially converts to photons (a_2).
#
# Validation (Boccaletti 1970, weak-field regime):
#   P(graviton -> photon) = sin^2(kappa * Bpeak * R * sqrt(pi/2))
# For defaults (kappa=1, Bpeak=0.1, R=5):
#   P ≈ sin^2(0.1 * 5 * sqrt(pi/2)) ≈ sin^2(0.627) ≈ 0.343
#
# Ref: Boccaletti et al. (1970, Nuovo Cimento 70B:129)
#      Domcke, Garcia-Cely & Lee (2025, arXiv:2507.16609)
#
# Running:
#   cd examples/gertsenshtein && bash run_localized.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/gertsenshtein_localized_output

echo "=== Gertsenshtein Effect: Localized B-field Scattering ==="
echo ""

# Step 1: Derive coupled equations from Einstein-Maxwell Lagrangian with Gaussian B(z)
echo "--- Step 1: Derive ---"
tidal derive theory_localized.toml
echo ""

# Step 2: Inspect the derived JSON (same 6 fields as uniform case; coefficients now z-dependent)
echo "--- Step 2: Inspect ---"
tidal inspect ../data/gertsenshtein_localized.json
echo ""

# Step 3: Simulate — graviton wavepacket scattering through localized B-field region
#
# Domain: [-100, 100] (200 units), 1024 grid points (dx ≈ 0.195)
# Non-periodic (Neumann BCs) — periodic BCs are incorrect here since Abar_y is not periodic
# Gaussian wavepacket centered at z = -50, width σ = 5, k = 2.0 (massless dispersion limit)
# t_end = 120: wavepacket travels 50 units to B-field center (t≈50), then 30 more units
#   out of the B-field region before hitting the right boundary (t≈150).
#   Stopping at 120 ensures exactly one passage through the B-field.
#
# Note: --what energy is NOT reliable with position-dependent coefficients (known limitation).
# Use --what conversion as the primary validation metric.
echo "--- Step 3: Simulate ---"
tidal simulate ../data/gertsenshtein_localized.json \
  --grid-shape 1024 \
  "--bounds=-100:100" \
  --ic gaussian \
  --ic-wavevector 2.0 \
  --ic-amplitude 0.1 \
  --ic-width 5.0 \
  "--ic-center=-50.0" \
  --ic-component h_7 \
  --t-end 120.0 \
  --param kappa=1.0 --param Bpeak=0.1 --param R=5.0 \
  --output "$OUT"

echo ""

# Step 4: Measure conversion probability
# Expected: P = sin^2(kappa * Bpeak * R * sqrt(pi/2))
#           = sin^2(1.0 * 0.1 * 5.0 * 1.2533) = sin^2(0.6267) ≈ 0.343
echo "--- Step 4: Measure ---"
tidal measure "$OUT" --what conversion --source h_7 --target a_2 \
  --param kappa=1.0 --param Bpeak=0.1 --param R=5.0
echo ""

# Step 5: Plots
echo "--- Step 5: Plots ---"
tidal plot "$OUT" --type heatmap --field h_7 \
  --title "Graviton h_+ scattering through Gaussian B(z)" \
  --output "$OUT/heatmap_h_7.png" --quiet

tidal plot "$OUT" --type heatmap --field a_2 \
  --title "Photon a_y (Gertsenshtein conversion, localized B)" \
  --output "$OUT/heatmap_a_2.png" --quiet

tidal plot "$OUT" --type amplitude \
  --fields h_7,a_2 \
  --title "Gertsenshtein h_+ → a_y (Gaussian B-field)" \
  --output "$OUT/amplitude.png" --quiet

echo ""
echo "=== Done ==="
echo "Plots saved to: $OUT/"
ls "$OUT"/*.png 2>/dev/null | sed 's|.*/|  |'
