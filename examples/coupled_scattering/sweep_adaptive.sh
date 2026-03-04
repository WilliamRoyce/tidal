#!/usr/bin/env bash
# Adaptive Coupling Sweep — Automatically Refines Near Sharp Features
#
# Demonstrates TIDAL's adaptive parameter sampling (F2a).  Starts with a
# coarse 5-point grid over coupling strength g0 and iteratively bisects
# intervals with high curvature in P_max until the budget is exhausted
# or all intervals are smooth.
#
# Expected behavior: points cluster near the instability knee (g0 ~ 0.7-0.9)
# where P_max transitions from weak to strong conversion.  Flat regions
# (g0 < 0.1) receive fewer points.
#
# Running this script:
#   cd examples/coupled_scattering && bash sweep_adaptive.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_adaptive

echo "=== Adaptive Coupling Sweep: g0 = 0.01 to 0.9, budget = 20 ==="

tidal sweep ../data/coupled_scattering.json \
  --sweep "g0=0.01:0.9:5" \
  --adaptive-metric P_max \
  --adaptive-budget 20 \
  --adaptive-threshold 0.01 \
  --param mPhi2=1.0 --param mChi2=1.0 --param R=5.0 \
  --measure conversion,mixing,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 128 \
  --bounds=-40:40 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 40.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo "  Note: adaptive sampling concentrates points near sharp features."
echo "  Compare with a uniform 20-point sweep to see the difference."
echo ""
echo "Plot the results:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
echo "  tidal plot $OUT/ --type sweep --metric L_mix"
