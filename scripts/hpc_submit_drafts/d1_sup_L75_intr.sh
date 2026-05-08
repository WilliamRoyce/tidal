#!/usr/bin/env bash
# Phase 6.J — D1 sup at L=75 INTR (bounds-dependence bracketing point, #344)
#
# Companion to d1_amp_L75_intr.sh. Fire ONLY after the amp job has
# completed (INTR has MaxSubmitPU=1).
#
# Phase 6.G observed a 3.4-nat shift in sup log Z between L=50 and L=100
# (+16.329 → +12.957, in opposite direction to amp). This run adds the
# L=75 bracketing point to verify monotonicity.
#
# Settings: matched to Phase 6.G sup midres (28967862, grid=128) and
# Phase 6.H sup intr-xreduced (28985713, grid=128/nlive=300/L=100).
# ic_wavevector=2π/75 (fundamental mode).
#
# Expected wall: ~30-50 min (sup landscape is multi-modal — 28985713 took
# 40 min at L=100; expect similar at L=75).
#
# Predicted log Z (sup): between +16.329 (L=50) and +12.957 (L=100), so
# ~+14.6 if linear in L. Anomalous if outside [+12, +17].
#
# Reference: docs/PHASE_6_COMPARISON.md §"Phase 6.G", #344

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_sup_L75_intr --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=uniform:-1:1" --prior "alpha2=uniform:-2:2" \
    --prior "alpha3=log_uniform:0.05:2" --prior "delta1=uniform:-2:2" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 128 --bounds=0:75 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.08377580409572781 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 300 \
    --num-repeats 5 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_sup_L75_intr'
