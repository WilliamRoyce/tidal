#!/usr/bin/env bash
# Phase 6.K — Stage A v5 sup at τ=0.15 (post-#341 probe)
#
# Phase 6.E reran the published Stage A v5 sup (28477675) at τ=0.3 (28859477)
# and got log Z = +0.602 ± 0.052 (Δ = −0.052 vs +0.654 original).
# This re-run repeats with τ=0.15 to close the probe-tightening cross-check.
#
# Settings: identical to 28859477 / 28477675. Single INTR, ~5 min wall.
# Reference: docs/PHASE_6_COMPARISON.md §"Phase 6.E"

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name dp_suppress_v5_tau015 --ntasks 76 \
  --cmd 'tidal sample examples/data/dark_photon_plasma.json \
    --prior "mA2=log_uniform:0.001:1.0" \
    --prior "deltam=uniform:-0.5:0.5" \
    --prior "xi=log_uniform:0.05:20.0" \
    --prior "alpha3=log_uniform:0.001:0.5" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 400 \
    --num-repeats 5 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/dp_suppress_v5_tau015'
