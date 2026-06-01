#!/usr/bin/env bash
# v3 Stage A sup — Dark-Photon-Plasma, tachyon-permissive, compactified priors
#
# Reference: 28477675 (Stage A v5 sup, log Z=+0.654±0.056); Phase 6.K rerun
# at τ=0.15 gave +0.715±0.062 (consistent — null robust under v2 probe).
# v3 re-run drops the probe gate and Hwang-Noh entirely.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name stage_a_sup_v3 --ntasks 76 \
  --cmd 'tidal sample examples/data/dark_photon_plasma.json \
    --prior "mA2=log_uniform:1e-3:1e3" \
    --prior "deltam=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/stage_a_sup_v3'
