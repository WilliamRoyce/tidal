#!/usr/bin/env bash
# Gertsenshtein with Massive Photon — 1D Resonance Scan (fast validation)
#
# Fixes B0 = 0.10 and scans mA² across the Raffelt-Stodolsky Lorentzian.
#
# Physics (Raffelt & Stodolsky 1988):
#   On resonance (mA² = κ²B₀²/2 ≈ 0.005): P_max → 1, amplitude
#   Off-resonance: P_max = (2·κ·B₀·k)² / [(2·κ·B₀·k)² + (mA² - κ²B₀²/2)²]
#   Lorentzian HWHM in mA² is 2·κ·B₀·k ≈ 0.4 for κ=1, B₀=0.10, k=2.0106.
#
# To see the suppression, the grid must extend well beyond coupling²:
#   mA² ∈ [0, 1.0] with 40 points (Δ=0.025)
#   At mA²=0.005 (on-resonance):   P_max ≈ 1.0
#   At mA²=0.405 (Δm²=coupling):   P_max ≈ 0.5
#   At mA²=0.805 (Δm²=2·coupling): P_max ≈ 0.2
#   This fully resolves the Lorentzian shape.
#
# Ref: Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#      Domcke, Garcia-Cely & Lee (2025), arXiv:2507.16609
#
# Running:
#   cd examples/gertsenshtein_proca && bash sweep_resonance_1d.sh

set -euo pipefail
cd "$(dirname "$0")"

JSON=../data/gertsenshtein_proca.json
OUT=../data/sweep_proca_resonance_1d

echo "=== Massive Photon Gertsenshtein: 1D Resonance Scan ==="
echo "    B0 = 0.10 fixed, mA2 in [0, 1.0] with 40 points"
echo "    Peak at mA2 = kappa^2 * B0^2 / 2 = 0.005, FWHM ~ 2*coupling ~ 0.8"
echo ""

tidal sweep "$JSON" \
  --sweep "mA2=0.0:1.0:40" \
  --param kappa=1.0 --param B0=0.10 \
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
echo "  tidal plot $OUT --type sweep --metric P_max \\"
echo "    --title 'Gertsenshtein P_max(mA2) — 1D resonance scan at B0=0.10' \\"
echo "    --overlay '4*(2.0106*1.0*0.10)^2/((mA2 - 0.005)^2 + 4*(2.0106*1.0*0.10)^2)'"
