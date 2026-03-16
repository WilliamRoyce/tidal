#!/usr/bin/env bash
# Coupled Scalars (Gertsenshtein Effective) — Full derive → inspect → simulate → measure pipeline
#
# Physics: Effective 1+1D graviton-photon mixing from Einstein-Maxwell.
# The gradient coupling B0*h*∂_x(a) transfers energy between the graviton
# mode (h) and photon mode (a), producing Rabi oscillations in conversion
# probability P(t) = sin²(κ*B0*t/2) for massless vacuum.
#
# Parameters:
#   kappa  = graviton-photon coupling (= √(16πG) in natural units)
#   B0     = background magnetic field strength
#   omegaP2 = plasma frequency squared (photon effective mass)
#   mg2    = graviton mass squared
#
# Refs:
#   Gertsenshtein (1962), JETP 14:84
#   Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#
# Running:
#   cd examples/coupled_scalars && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system (2 fields: h_0, a_0)
tidal inspect ../data/coupled_scalars.json

# ============================================================================
# Run 1: Graviton wavepacket → photon conversion (massless vacuum)
# ============================================================================
# Plane-wave IC in h_0 only; gradient coupling converts h → a.
# With B0=0.1, kappa=1: conversion period T = 2π/(κB₀) ≈ 63, t_end=50 covers ~80%.
tidal simulate ../data/coupled_scalars.json \
  --param kappa=1.0 --param B0=0.1 --param omegaP2=0.0 --param mg2=0.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-component h_0 \
  --ic-wavevector 2.0 \
  --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output ../data/coupled_scalars_output

# --- Run 1 measurements ---
# P(t) = E_a(t) / E_h(0): what fraction of graviton energy converted to photon
tidal measure ../data/coupled_scalars_output \
  --what conversion,conservation \
  --source h_0 --target a_0 \
  --param kappa=1.0 --param B0=0.1 --param omegaP2=0.0 --param mg2=0.0

# Dispersion relation omega(k)
tidal measure ../data/coupled_scalars_output \
  --what dispersion \
  --source h_0 \
  --param kappa=1.0 --param B0=0.1 --param omegaP2=0.0 --param mg2=0.0

# --- Run 1 plots ---
tidal plot ../data/coupled_scalars_output --type heatmap --field h_0 \
  --title "Graviton h (Gertsenshtein conversion)" --quiet
tidal plot ../data/coupled_scalars_output --type heatmap --field a_0 \
  --title "Photon a (converted from h)" --quiet
tidal plot ../data/coupled_scalars_output --type amplitude --quiet
tidal plot ../data/coupled_scalars_output --type energy --quiet
tidal plot ../data/coupled_scalars_output --type conservation --quiet

echo ""
echo "=== Run complete ==="
echo "  Output: ../data/coupled_scalars_output/"
echo "  Expected: P_max = sin²(κ*B0*t_peak/2) = sin²(π/2) ≈ 1.0 at t ≈ π/(κB₀) ≈ 31.4"
