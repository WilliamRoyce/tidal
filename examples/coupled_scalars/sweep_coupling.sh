#!/usr/bin/env bash
# Coupled Scalars — Sweep coupling strength gCpl
#
# Sweeps the coupling constant from weak to strong coupling and measures
# conversion probability P(t) = E_chi(t) / E_phi(0) at each point.
# Uses a standing Gaussian wave packet (no wavevector) in a homogeneous
# coupling field, so conversion occurs via uniform Rabi-like oscillations.
#
# Expected physics:
#   - P_max increases with gCpl (stronger coupling = more energy transfer)
#   - P(t) oscillates with mixing frequency proportional to gCpl
#   - Energy conservation should be tight (max |dE/E| < 1e-3)
#
# Running this script:
#   cd examples/coupled_scalars && bash sweep_coupling.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_sweep

echo "=== Coupled Scalars Coupling Sweep: gCpl = 0.01 to 1.0 (8 points, log scale) ==="

tidal sweep ../data/coupled_scalars.json \
  --sweep "gCpl=0.01:1.0:8:log" \
  --param mPhi2=1.0 --param mChi2=4.0 \
  --measure conversion,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center 15.0 --ic-width 3.0 \
  --t-end 10.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Plot the results:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
echo "  tidal plot $OUT/ --type sweep --metric P_max,max_energy_error"
