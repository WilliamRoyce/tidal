#!/usr/bin/env bash
# Gertsenshtein with Massive Photon — 2D Resonance Map
#
# Maps conversion probability over B0 x mA2, revealing the resonance
# structure. The graviton acquires an effective mass mg2_eff = kappa^2*B0^2/2
# from the background EM energy. Maximum conversion occurs when the
# graviton and photon effective masses match.
#
# Ref: Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#      Domcke, Garcia-Cely & Lee (2025), arXiv:2507.16609
#
# Running:
#   cd examples/gertsenshtein_proca && bash sweep_resonance.sh

set -euo pipefail
cd "$(dirname "$0")"

JSON=../data/gertsenshtein_proca.json
OUT=../data/sweep_proca_resonance

echo "=== Massive Photon Gertsenshtein: 2D Resonance Map ==="
echo "    B0 = 0.01 to 0.25 (10 pts) x mA2 = 0.0 to 2.0 (30 pts)"
echo "    300 simulations total"
echo ""
echo "Physics (Raffelt & Stodolsky 1988):"
echo "  Peak at mA² = κ²B₀²/2 ∈ [5e-5, 0.031], but the Lorentzian HWHM"
echo "  in mA² is 2·κ·B₀·k ~ 0.04-0.81 across the B0 range. The grid"
echo "  must extend past 2·(HWHM at B0=0.25) ~ 1.6 to resolve the fall-off."
echo ""
echo "Predicted resonance line: mA2 = kappa^2 * B0^2 / 2"
echo ""

tidal sweep "$JSON" \
  --sweep "B0=0.01:0.25:10" \
  --sweep "mA2=0.0:2.0:30" \
  --param kappa=1.0 \
  --measure conversion,conservation \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output "$OUT" --resume

echo ""
echo "Results: $OUT/"
echo ""
echo "Visualize:"
echo "  tidal plot $OUT/ --type sweep-scatter --x-param B0 --y-param mA2 --metric P_max"
