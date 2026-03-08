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
# After constraint elimination (traceless, transverse → gradient-zero), the
# 3+1D system reduces to 6 fields in 1+1D:
#   h_7 ↔ a_2  (plus polarization — Gertsenshtein channel)
#   h_5 ↔ a_1  (cross polarization — decoupled, stays zero)
#   a_0, a_3   (gauge modes — decoupled free waves, stay zero)
#
# Initial condition: monochromatic plane wave on h_7 (graviton h_+).
# Plane-wave IC matches the analytical derivation (Palessandro & Rothman
# 2023, Eq. 26) and avoids spectral broadening of a Gaussian wave packet.
# Wavevector k = 2π·32/100 ≈ 2.011 ensures k >> B₀κ (massless limit)
# and integer wavelengths in the domain (no spectral leakage).
#
# Validation: P(graviton -> photon) = sin^2(kappa * B0 * D / 2)
# with D = c * t_end = 200.
# Derived from eigenmode analysis of coupled h_7/a_2 system.
#
# Running:
#   cd examples/gertsenshtein && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/gertsenshtein_output

echo "=== Gertsenshtein Effect: Graviton-Photon Conversion ==="
echo ""

# Step 1: Derive coupled equations from Einstein-Maxwell Lagrangian
echo "--- Step 1: Derive ---"
tidal derive theory.toml
echo ""

# Step 2: Inspect the derived JSON (expect 6 fields after elimination)
echo "--- Step 2: Inspect ---"
tidal inspect ../data/gertsenshtein.json
echo ""

# Step 3: Simulate — plane wave graviton propagating through B0 region
# k = 2π·32/100 ≈ 2.011 (32 wavelengths in domain, 16 pts/wavelength)
# t_end = 200 → D = 200, giving ~2 full Rabi oscillation cycles at B0=0.3
echo "--- Step 3: Simulate ---"
tidal simulate ../data/gertsenshtein.json \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-wavevector 2.0106 \
  --ic-amplitude 0.1 \
  --ic-component h_7 \
  --t-end 200.0 \
  --param kappa=1.0 --param B0=0.3 \
  --output "$OUT"

echo ""

# Step 4: Measure energy conservation and conversion
echo "--- Step 4: Measure ---"
tidal measure "$OUT" --what conservation --param kappa=1.0 --param B0=0.3
tidal measure "$OUT" --what conversion --param kappa=1.0 --param B0=0.3
tidal measure "$OUT" --what energy --param kappa=1.0 --param B0=0.3
echo ""

# Step 5: Plots
echo "--- Step 5: Plots ---"
tidal plot "$OUT" --type heatmap --field h_7 \
  --title "Graviton h_+ (plus polarization)" \
  --output "$OUT/heatmap_h_7.png" --quiet

tidal plot "$OUT" --type heatmap --field a_2 \
  --title "Photon a_y (Gertsenshtein conversion)" \
  --output "$OUT/heatmap_a_2.png" --quiet

tidal plot "$OUT" --type amplitude \
  --fields h_7,a_2 \
  --title "Gertsenshtein h_+ → a_y conversion" \
  --output "$OUT/amplitude.png" --quiet

echo ""
echo "=== Done ==="
echo "Plots saved to: $OUT/"
ls "$OUT"/*.png 2>/dev/null | sed 's|.*/|  |'
