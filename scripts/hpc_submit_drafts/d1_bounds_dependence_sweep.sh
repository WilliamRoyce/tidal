#!/usr/bin/env bash
# Phase 6.J — Bounds-dependence sweep (#344) at τ=0.15
#
# Phase 6.G observed a 3-nat shift in both amp and sup log Z between
# L=50 and L=100 (in opposite directions: amp +2.8, sup −3.4). This
# sweep brackets L=100 with L=75 and L=150 to verify the trend is
# monotonic and L=100 is in a stable regime, not a local extremum.
#
# 4 paired chains at xreduced INTR fidelity (matches 28985713):
#   grid=128, nlive=300, snapshots=2, num-repeats=5
#   ic_wavevector = 2π/L (fundamental mode)
# Each ~30-40 min INTR; 4 jobs total.
#
# Reference: docs/PHASE_6_COMPARISON.md §"Phase 6.G", #344

set -euo pipefail

submit_pair () {
  local L=$1
  local k_ic
  k_ic=$(python3 -c "import math; print(2*math.pi/$L)")

  for mode in amp sup; do
    if [ "$mode" = "amp" ]; then likely="P_max:maximize"; else likely="P_max:minimize"; fi
    bash scripts/hpc_shuttle.sh submit \
      --template scripts/hpc_templates/polychord_intr.sbatch \
      --name "d1_${mode}_L${L}_intr" --ntasks 76 \
      --cmd "tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
        --prior \"alpha1=uniform:-1:1\" --prior \"alpha2=uniform:-2:2\" \
        --prior \"alpha3=log_uniform:0.05:2\" --prior \"delta1=uniform:-2:2\" \
        --likelihood \"$likely\" \
        --baseline-formula \"sin(kappa*B0*t_end/2)**2\" \
        --param kappa=1.0 --param B0=0.01 \
        --grid-shape 128 --bounds=0:$L --periodic \
        --ic plane-wave --ic-component h_5 --ic-wavevector $k_ic --ic-amplitude 1e-2 \
        --t-end 10 --snapshots 2 \
        --measure conversion,peak_conversion --source h_5 --target a_1 \
        --method nested --sampler polychord --nlive 300 \
        --num-repeats 5 --precision-criterion 0.01 \
        --output \${RESULTS_DIR}/d1_${mode}_L${L}_intr"
  done
}

submit_pair 75
submit_pair 150
