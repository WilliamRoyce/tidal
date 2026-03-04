#!/usr/bin/env bash
# 1D Coupling Strength Sweep — Gertsenshtein Conversion vs g0
#
# Sweeps the coupling constant g0 from 0.1 to 0.9 (10 points) with fixed
# masses and interaction length.  Measures conversion probability and mixing
# length at each point.
#
# Expected physics: P_max increases monotonically with g0 (stronger coupling
# = more energy transfer from phi to chi).  Near the stability threshold
# (g0^2 → mPhi2*mChi2 = 1), P_max approaches O(1).
#
# Running this script:
#   cd examples/coupled_scattering && bash sweep.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_coupling

echo "=== Coupling Strength Sweep: g0 = 0.1 to 0.9 (10 points) ==="

tidal sweep ../data/coupled_scattering.json \
  --sweep "g0=0.1:0.9:10" \
  --param mPhi2=1.0 --param mChi2=1.0 --param R=5.0 \
  --measure conversion,mixing,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 256 \
  --bounds=-40:40 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 60.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo "  sweep.json     — sweep definition + provenance"
echo "  results.csv    — self-contained metrics table"
echo "  results.json   — same data in JSON format"
echo ""
echo "Plot the results:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
echo "  tidal plot $OUT/ --type sweep --metric L_mix"
echo "  tidal plot $OUT/ --type sweep-compare --metric conversion"
