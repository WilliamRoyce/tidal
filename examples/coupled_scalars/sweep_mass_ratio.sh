#!/usr/bin/env bash
# Mass Ratio Sweep — Resonance Structure of Coupled Scalars
#
# Sweeps the target mass mChi2 from 0.1 to 10.0 at fixed coupling,
# demonstrating the resonance peak when mChi2 ≈ mPhi2 = 1.0
# (degenerate masses = exact frequency matching at all k).
#
# Measures conversion probability, velocity, and resonance analysis
# to show how conversion bandwidth narrows as the mass ratio moves
# away from resonance.
#
# Expected physics:
#   - P_max peaks sharply at mChi2 ≈ 1.0 (resonant)
#   - Conversion bandwidth is widest at resonance, narrowing for |mChi2 - mPhi2| >> 0
#   - Group velocity mismatch increases off-resonance (decoherence)
#   - v_group = k/sqrt(k^2 + m^2) differs for fields with different masses
#
# Uses --ic-wavevector 3 to create a propagating wave packet, which
# concentrates Fourier energy around k=3 for cleaner velocity/resonance
# analysis.
#
# Demonstrates features F3 (velocity/resonance) and F6 (spectrum scalars).
#
# Running this script:
#   cd examples/coupled_scalars && bash sweep_mass_ratio.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_mass_ratio

echo "=== Mass Ratio Sweep: mChi2 = 0.1 to 10.0 (15 points, log scale) ==="

tidal sweep ../data/coupled_scalars.json \
  --sweep "mChi2=0.1:10.0:15:log" \
  --param mPhi2=1.0 --param gCpl=0.3 \
  --measure conversion,velocity,resonance,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 128 \
  --bounds 0:80 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center 20.0 --ic-width 4.0 --ic-wavevector 3 \
  --t-end 20.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Plot resonance structure:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
echo "  tidal plot $OUT/ --type sweep --metric conversion_bandwidth"
echo "  tidal plot $OUT/ --type sweep --metric n_resonant_modes"
