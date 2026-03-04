#!/usr/bin/env bash
# Grid Convergence Study — Energy Conservation vs Resolution
#
# Runs the same physics (resonant conversion, g0=0.5) at grid sizes
# 32, 64, 128, 256, 512.  Measures energy conservation at each
# resolution to determine the spatial convergence order.
#
# Expected physics: max |dE/E| decreases as O(h^2) for leapfrog with
# 2nd-order spatial operators.  The convergence plot should show a
# straight line with slope ~2 on a log-log scale.
#
# Running this script:
#   cd examples/coupled_scattering && bash convergence.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_convergence

echo "=== Grid Convergence Study: N = 32, 64, 128, 256, 512 ==="

tidal sweep ../data/coupled_scattering.json \
  --converge "32,64,128,256,512" \
  --param mPhi2=1.0 --param mChi2=1.0 --param g0=0.5 --param R=5.0 \
  --measure conservation \
  --bounds=-40:40 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 40.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Plot the convergence:"
echo "  tidal plot $OUT/ --type convergence --metric max_energy_error"
