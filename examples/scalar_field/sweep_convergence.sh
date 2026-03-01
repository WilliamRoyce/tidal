#!/usr/bin/env bash
# Grid Convergence Study — Klein-Gordon Scalar Field
#
# Runs a convergence study on the simplest TIDAL example (free scalar field)
# to verify spatial convergence order.  Measures energy conservation at
# grid sizes 32, 64, 128, 256.
#
# Demonstrates that the sweep framework works for any theory, including
# single-field systems without cross-field coupling.
#
# Expected physics: max |dE/E| decreases as O(h^2) for leapfrog with
# second-order spatial operators.
#
# Running this script:
#   cd examples/scalar_field && bash sweep_convergence.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_convergence_kg

echo "=== Klein-Gordon Convergence Study: N = 32, 64, 128, 256 ==="

tidal sweep ../data/klein_gordon_1d.json \
  --converge "32,64,128,256" \
  --param m2=1.0 \
  --measure conservation \
  --bounds 0:50 \
  --periodic \
  --ic gaussian --ic-center 25.0 --ic-width 3.0 \
  --t-end 20.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Plot the convergence:"
echo "  tidal plot $OUT/ --type convergence --metric max_energy_error"
