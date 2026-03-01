#!/usr/bin/env bash
# 2D Parameter Sweep — Coupling Strength x Mass Ratio
#
# Sweeps g0 (coupling) and mChi2 (target mass) simultaneously, creating
# a 2D grid of simulations.  Produces a heatmap of conversion probability
# over the (g0, mChi2) parameter plane.
#
# Expected physics: conversion is strongest on the resonant diagonal
# (mChi2 ≈ mPhi2 = 1.0) and increases with coupling strength.  Off-resonant
# masses suppress conversion due to momentum mismatch.
#
# Running this script:
#   cd examples/coupled_scattering && bash sweep_2d.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_coupling_mass_2d

echo "=== 2D Sweep: g0 x mChi2 (5 x 5 = 25 runs) ==="

tidal sweep ../data/coupled_scattering.json \
  --sweep "g0=0.1:0.9:5" \
  --sweep "mChi2=0.5:4.0:5" \
  --param mPhi2=1.0 --param R=5.0 \
  --measure conversion \
  --source phi_0 --target chi_0 \
  --grid-shape 128 \
  --bounds=-40:40 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 40.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Plot the 2D heatmap:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
