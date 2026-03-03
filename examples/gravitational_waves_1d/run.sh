#!/usr/bin/env bash
# Gravitational Waves (TT gauge + z-propagation) — Reduced 1+1D pipeline
#
# Physics: Linearized gravity in TT gauge with plane-wave reduction along z.
# The full 3+1D linearized Einstein equations (10 symmetric h_ab components,
# reduced to 5 after temporal+traceless+transverse gauge conditions) further
# reduce to just 2 uncoupled 1D wave equations for the plus and cross
# polarizations:
#   d2t(h_+) = d2z(h_+)     (plus polarization, h_xx = -h_yy)
#   d2t(h_x) = d2z(h_x)     (cross polarization, h_xy)
#
# This is the primary use case for plane-wave reduction: testing GW
# propagation along a single spatial axis at dramatically lower cost
# than a full 3D simulation.
#
# Note: The derived JSON contains 12 fields (5 dynamical + 7 constraints).
# IDA correctly handles the TT gauge constraints. The physical DOFs are
# h_4 (h_+ = h_xx) and h_5 (h_x = h_xy), each satisfying d2t = d2z.
#
# The unified constraint IC solver automatically determines:
#   h_7 = -h_4 (from transverse_z constraint: gradient(h_4) + gradient(h_7) = 0)
#   h_9 = -(h_4 + h_7) (from traceless constraint)
#   h_6 = h_8 = 0 (from transverse_x/y constraints)
#
# Running this script:
#   cd examples/gravitational_waves_1d && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/gw_plane_wave_1d_output

# Derive equations from Lagrangian (requires wolframscript)
# tidal derive theory.toml

# Inspect the equation system (should show 1+1D, 12 fields with constraints)
tidal inspect ../data/gw_plane_wave_1d.json

# Run 1D simulation (Gaussian pulse on h_+ polarization, periodic BCs)
tidal simulate ../data/gw_plane_wave_1d.json \
  --grid-shape 256 \
  --bounds=-10:10 \
  --periodic \
  --ic gaussian \
  --ic-width 1.0 \
  --ic-component h_4 \
  --t-end 8.0 \
  --output "$OUT"

# --- Analysis plots ---

# Spacetime heatmap: h_4 (plus polarization) evolution
tidal plot "$OUT" --type heatmap --field h_4 \
  --title "GW h_+ polarization (x-t)" \
  --output "$OUT/heatmap_h_4.png" --quiet

# Spacetime heatmap: h_7 (traceless partner, should be -h_4)
tidal plot "$OUT" --type heatmap --field h_7 \
  --title "GW h_7 = -h_+ (traceless partner)" \
  --output "$OUT/heatmap_h_7.png" --quiet

# Profile evolution: h_4 at multiple time snapshots
tidal plot "$OUT" --type profile --field h_4 \
  --time-indices 0,25,50,75,100 \
  --title "h_+ profile evolution" \
  --output "$OUT/profile_h_4.png" --quiet

# Profile evolution: h_7 at multiple time snapshots
tidal plot "$OUT" --type profile --field h_7 \
  --time-indices 0,25,50,75,100 \
  --title "h_7 profile evolution" \
  --output "$OUT/profile_h_7.png" --quiet

# Multi-field peak amplitude
tidal plot "$OUT" --type amplitude \
  --fields h_4,h_5,h_7 \
  --title "GW field amplitudes" \
  --output "$OUT/amplitude_dynamical.png" --quiet

# Compare h_4 initial/final vs h_7
tidal plot "$OUT" --type compare \
  --fields h_4,h_7 \
  --title "h_4 vs h_7 (TT constraint: h_7 = -h_4)" \
  --output "$OUT/compare_h4_h7.png" --quiet

# Energy conservation measurement (currently blocked by unsupported
# mixed_* Hamiltonian operators — see .github-issues-pending.md)
# tidal plot "$OUT" --type conservation
# tidal plot "$OUT" --type hamiltonian --fields h_4,h_5,h_7

echo ""
echo "=== Run 1 (static Gaussian) complete. Plots: $OUT/ ==="

# ============================================================================
# Run 2: Unidirectional GW pulse (right-mover via formula + velocity)
# ============================================================================
# The static Gaussian above splits into left+right movers. Using
# --ic-formula-velocity, we lock field and velocity together to create
# a clean unidirectional GW burst — the physically correct way to
# initialize a gravitational wave pulse.
#
# For a right-mover: h_+ = envelope*cos(kx), dh_+/dt = +k*envelope*sin(kx)
# The constraint solver automatically sets h_7 = -h_4 (traceless partner).

OUT2=../data/gw_plane_wave_1d_travelling

tidal simulate ../data/gw_plane_wave_1d.json \
  --grid-shape 256 \
  --bounds=-10:10 \
  --periodic \
  --ic formula \
  --ic-formula 'exp(-x**2 / 2) * cos(3 * x)' \
  --ic-formula-velocity '3 * exp(-x**2 / 2) * sin(3 * x)' \
  --ic-component h_4 \
  --t-end 8.0 \
  --output "$OUT2"

# Heatmap: clean rightward propagation (compare vs symmetric Run 1)
tidal plot "$OUT2" --type heatmap --field h_4 \
  --title "GW h_+ travelling pulse (right-mover)" \
  --output "$OUT2/heatmap_h_4.png" --quiet

tidal plot "$OUT2" --type heatmap --field h_7 \
  --title "GW h_7 = -h_+ travelling (traceless partner)" \
  --output "$OUT2/heatmap_h_7.png" --quiet

tidal plot "$OUT2" --type profile --field h_4 \
  --time-indices 0,25,50,75,100 \
  --title "h_+ travelling profile evolution" \
  --output "$OUT2/profile_h_4.png" --quiet

tidal plot "$OUT2" --type amplitude \
  --fields h_4,h_5,h_7 \
  --title "GW amplitudes (travelling pulse)" \
  --output "$OUT2/amplitude.png" --quiet

echo ""
echo "=== Run 2 (travelling pulse) complete ==="
echo ""
echo "All runs complete. Generated plots:"
echo ""
echo "  Static Gaussian ($OUT/):"
ls "$OUT"/*.png 2>/dev/null | sed 's|.*/|    |'
echo ""
echo "  Travelling pulse ($OUT2/):"
ls "$OUT2"/*.png 2>/dev/null | sed 's|.*/|    |'
