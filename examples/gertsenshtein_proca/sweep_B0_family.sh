#!/usr/bin/env bash
# Gertsenshtein with Massive Photon — B0 Family Sweep
#
# Sweeps the background magnetic field strength B0 at multiple
# photon mass values mA2, producing a family of P(B0) curves
# that show how mass detuning suppresses and shifts the Rabi oscillation.
#
# Expected physics:
#   mA2=0: P(B0) = sin²(κ*B0*t_end/2) — standard Gertsenshtein
#   mA2>0: suppressed and shifted oscillation via Raffelt-Stodolsky
#
# Ref: Gertsenshtein (1962), JETP 14:84
#      Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#
# Running:
#   cd examples/gertsenshtein_proca && bash sweep_B0_family.sh

set -euo pipefail
cd "$(dirname "$0")"

JSON=../data/gertsenshtein_proca.json

echo "=== Massive Photon Gertsenshtein: B0 Family Sweep ==="
echo "    B0 = 0.005 to 0.25 (30 points), multiple mA2 values"
echo ""

for mA2 in 0.0 0.2 0.5 1.0 2.0; do
  OUT=../data/sweep_proca_B0_mA2${mA2}

  echo "--- mA2 = $mA2 ---"

  tidal sweep "$JSON" \
    --sweep "B0=0.005:0.25:30" \
    --param kappa=1.0 --param mA2="$mA2" \
    --measure conversion,conservation \
    --source h_7 --target a_2 \
    --grid-shape 512 \
    --bounds 0:100 \
    --periodic \
    --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 --ic-amplitude 0.1 \
    --t-end 50.0 \
    --output "$OUT" --resume

  # For mA2=0, overlay the analytical Gertsenshtein formula
  if [ "$mA2" = "0.0" ]; then
    tidal plot "$OUT" --type sweep \
      --metric P_final \
      --title "P(B0) — mA2=$mA2 (massless)" \
      --overlay 'sin(kappa * B0 * t_end / 2)**2' \
      --output "$OUT/sweep.png" --quiet
  else
    tidal plot "$OUT" --type sweep \
      --metric P_max \
      --title "P(B0) — mA2=$mA2" \
      --output "$OUT/sweep.png" --quiet
  fi

  echo "  Results: $OUT/"
  echo ""
done

echo "Done. Visualize with:"
echo "  tidal plot ../data/sweep_proca_B0_*/ --type sweep --metric P_max"
