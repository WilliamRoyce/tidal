#!/usr/bin/env bash
# Gertsenshtein with Massive Photon — Detuning Sweep
#
# Sweeps the photon effective mass mA2 (= plasma frequency squared)
# at multiple B0 values, demonstrating detuning suppression of conversion.
#
# Expected physics (Raffelt-Stodolsky master formula):
#   P = sin²(2θ) × sin²(Δ_osc * D / 2)
#   tan(2θ) = κ*B0 / |Δ|,  Δ = -mA2/(2ω),  Δ_osc = √(Δ² + (κB0)²)
#
# At mA2=0 (massless): full Rabi oscillation P = sin²(κB0D/2)
# At mA2>>0: suppressed as P ~ (κB0*ω/mA2)² — detuned
#
# Ref: Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#      Domcke, Garcia-Cely & Lee (2025), arXiv:2507.16609, Sec. 6
#
# Running:
#   cd examples/gertsenshtein_proca && bash sweep_detuning.sh

set -euo pipefail
cd "$(dirname "$0")"

JSON=../data/gertsenshtein_proca.json

echo "=== Massive Photon Gertsenshtein: Detuning Sweep ==="
echo "    mA2 = 0.0 to 5.0 (40 points), multiple B0 values"
echo ""

for B0 in 0.05 0.1 0.15 0.2; do
  OUT=../data/sweep_proca_detuning_B0${B0}

  echo "--- B0 = $B0 ---"

  tidal sweep "$JSON" \
    --sweep "mA2=0.0:5.0:40" \
    --param kappa=1.0 --param B0="$B0" \
    --measure conversion,conservation \
    --source h_7 --target a_2 \
    --grid-shape 512 \
    --bounds 0:100 \
    --periodic \
    --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 --ic-amplitude 0.1 \
    --t-end 50.0 \
    --output "$OUT" --resume

  tidal plot "$OUT" --type sweep \
    --metric P_max \
    --title "Detuning P(mA2) — B0=$B0" \
    --output "$OUT/sweep.png" --quiet

  echo "  Results: $OUT/"
  echo ""
done

echo "Done. Visualize with:"
echo "  tidal plot ../data/sweep_proca_detuning_B0*/ --type sweep --metric P_max"
