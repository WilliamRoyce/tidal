#!/usr/bin/env bash
# Cylindrical Klein-Gordon (axial plane wave) — Reduced 1+1D pipeline
#
# Physics: Klein-Gordon in cylindrical coordinates, plane-wave reduction
# along z-axis (axial propagation). Since the z-axis in cylindrical
# coordinates is flat (g_zz = 1, no Christoffels), this reduces to
# a simple 1D massive wave equation: d2t(phi) = d2z(phi) - m2*phi.
# No volume element (flat axis), no position-dependent coefficients.
#
# This tests that reduction correctly eliminates all curved-coordinate
# complexity when the propagation axis is flat.
#
# Running this script:
#   cd examples/cylindrical_kg_1d && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/cylindrical_kg_1d_output

# Derive equations from Lagrangian (requires wolframscript)
# tidal derive theory.toml

# Inspect the equation system (should show 1+1D, 1 field)
tidal inspect ../data/cylindrical_kg_1d.json

# Run 1D simulation (Gaussian pulse on flat z-axis)
tidal simulate ../data/cylindrical_kg_1d.json \
  --grid-shape 128 \
  --bounds=-5:5 \
  --periodic \
  --ic gaussian \
  --ic-width 0.5 \
  --t-end 10.0 \
  --output "$OUT"

# --- Analysis plots ---

# Spacetime heatmap: phi_0 evolution (dispersive wave with mass)
tidal plot "$OUT" --type heatmap --field phi_0 \
  --title "Cylindrical KG (z-propagation, x-t)" \
  --output "$OUT/heatmap_phi_0.png"

# Profile evolution at multiple times
tidal plot "$OUT" --type profile --field phi_0 \
  --time-indices 0,25,50,75,100 \
  --title "phi_0 profile evolution" \
  --output "$OUT/profile_phi_0.png"

# Peak amplitude vs time (should show dispersive decay)
tidal plot "$OUT" --type amplitude --fields phi_0 \
  --title "Peak amplitude (massive dispersion)" \
  --output "$OUT/amplitude_phi_0.png"

# Hamiltonian energy decomposition
tidal plot "$OUT" --type hamiltonian --fields phi_0 \
  --title "Hamiltonian energy" \
  --output "$OUT/hamiltonian_phi_0.png"

# Energy conservation check
tidal plot "$OUT" --type conservation \
  --title "Energy conservation" \
  --output "$OUT/conservation.png"

# Compare initial vs final profile
tidal plot "$OUT" --type compare --fields phi_0 \
  --title "Initial vs final profile" \
  --output "$OUT/compare_phi_0.png"

echo ""
echo "Analysis complete. Plots saved to: $OUT/"
ls -la "$OUT"/*.png 2>/dev/null || echo "  (no plots generated)"
