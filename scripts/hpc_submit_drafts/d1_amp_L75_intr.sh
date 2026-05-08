#!/usr/bin/env bash
# Phase 6.J — D1 amp at L=75 INTR (bounds-dependence bracketing point, #344)
#
# Phase 6.G observed a 3-nat shift in amp log Z between L=50 and L=100
# (−2.091 → +0.679, +2.8 nats). This run adds a single bracketing point
# at L=75 to verify monotonicity and that L=100 isn't a local extremum.
#
# Settings: matched to Phase 6.G amp midres (28976470, grid=128) and
# Phase 6.H sup intr-xreduced (28985713, grid=128/nlive=300) for direct
# comparison. ic_wavevector=2π/75 (fundamental mode).
#
# Pairs with: d1_sup_L75_intr.sh — fire that AFTER this job completes
# (INTR has MaxSubmitPU=1 — only one INTR job per user at a time).
#
# Expected wall: ~30-40 min (bracketed by 28976470's 4 min at L=50/grid=128
# amp single-cluster, and 28985713's 40 min at L=100/grid=128 sup multi-cluster).
#
# Predicted log Z (amp): between −2.091 (L=50) and +0.679 (L=100), so ~−0.7
# if the trend is roughly linear in L. Anomalous if outside [−2.5, +1.0].
#
# Reference: docs/PHASE_6_COMPARISON.md §"Phase 6.G", #344

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_amp_L75_intr --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=uniform:-1:1" --prior "alpha2=uniform:-2:2" \
    --prior "alpha3=log_uniform:0.05:2" --prior "delta1=uniform:-2:2" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 128 --bounds=0:75 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.08377580409572781 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 300 \
    --num-repeats 5 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_amp_L75_intr'
