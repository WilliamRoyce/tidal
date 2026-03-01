#!/usr/bin/env bash
# Scattering Observables Sweep — Transmission & Reflection vs Coupling
#
# Sweeps coupling strength g0 and measures asymptotic scattering observables:
# P_transmitted (forward-propagating conversion) and P_reflected (backward
# scattering).  Uses the asymptotic measurement which decomposes the target
# field energy into directional components based on spectral flux.
#
# Expected physics:
#   - P_transmitted increases with g0 (more forward conversion)
#   - P_reflected remains small for smooth coupling profiles (Gaussian G(z))
#   - P_transmitted + P_reflected = P_final (total asymptotic conversion)
#
# Running this script:
#   cd examples/coupled_scattering && bash sweep_scattering.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_scattering

echo "=== Scattering Observables Sweep: g0 = 0.05 to 0.8 (12 points) ==="

tidal sweep ../data/coupled_scattering.json \
  --sweep "g0=0.05:0.8:12" \
  --param mPhi2=1.0 --param mChi2=1.0 --param R=5.0 \
  --measure asymptotic,peak_conversion,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 256 \
  --bounds=-40:40 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center=-20.0 --ic-width 3.0 --ic-wavevector 3 \
  --t-end 60.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Plot scattering observables:"
echo "  tidal plot $OUT/ --type sweep --metric P_transmitted"
echo "  tidal plot $OUT/ --type sweep --metric P_max,P_transmitted,P_reflected"
