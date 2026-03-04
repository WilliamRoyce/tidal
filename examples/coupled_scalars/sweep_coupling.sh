#!/usr/bin/env bash
# Coupled Scalars — Comprehensive coupling strength investigation
#
# Sweeps the coupling constant gCpl from weak to strong coupling,
# measuring conversion probability P(t), mixing timescale, asymptotic
# efficiency, effective mass renormalization, and energy conservation.
#
# Uses a standing Gaussian wave packet (no wavevector) in a homogeneous
# coupling field. The k-mode spread in the wave packet causes decoherence
# that averages Rabi oscillations into smooth sigmoid saturation curves.
#
# Expected physics:
#   - Weak coupling (gCpl << |m1² - m2²|/2):
#       P_max ≈ 4g²/(Δm²)² (perturbative), slow saturation
#   - Moderate coupling (gCpl ~ |m1² - m2²|/2 = 1.5):
#       P_max approaches 1, fastest energy transfer
#   - Strong coupling (gCpl >> |m1² - m2²|/2):
#       Normal mode splitting modifies dispersion; group velocity
#       mismatch between k-modes reduces coherent conversion.
#       Non-monotonic P_max may appear due to over-coupling.
#   - Rabi mixing angle: θ = ½ atan(2g / |m1² - m2²|)
#   - Homogeneous limit: P_max = sin²(2θ) = 4g²/[(Δm²)² + 4g²]
#   - Mixing length L_mix decreases with coupling (faster oscillations)
#   - Energy conservation tight: max |dE/E| < 1e-3
#
# Running this script:
#   cd examples/coupled_scalars && bash sweep_coupling.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_sweep

echo "=== Coupled Scalars Coupling Sweep: gCpl = 0.05 to 3.0 (20 points, log scale) ==="
echo "    mPhi2 = 1.0, mChi2 = 4.0  (mass gap Δm² = 3.0)"
echo "    Standing Gaussian on phi_0, domain [0,80], N=128, t_end=60"
echo ""

tidal sweep ../data/coupled_scalars.json \
  --sweep "gCpl=0.05:3.0:20:log" \
  --param mPhi2=1.0 --param mChi2=4.0 \
  --measure conversion,peak_conversion,mixing,asymptotic,effective_mass,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 128 \
  --bounds 0:80 \
  --periodic \
  --ic gaussian --ic-component phi_0 --ic-center 20.0 --ic-width 3.0 \
  --t-end 60.0 \
  --output "$OUT"

echo ""
echo "Results written to: $OUT/"
echo ""
echo "Visualize the results:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
echo "  tidal plot $OUT/ --type sweep --metric P_max,mixing_length"
echo "  tidal plot $OUT/ --type sweep-compare --measurement conversion"
echo "  tidal plot $OUT/ --type sweep --metric P_final,P_transmitted,P_reflected"
echo "  tidal plot $OUT/ --type sweep --metric m2_eff"
