#!/usr/bin/env bash
# B.4-full-sup — Stage A v3 sup PUBLICATION (standard queue, hi-res)
#
# Standard-queue version of stage_a_sup.sh. P_max:minimize for the
# suppression chain. Larger grid + nlive + num-repeats for
# publication-grade convergence. v3.1 priors (see B.4-full-amp).
#
# Reference: B.4b INTR smoke (29199129, log Z=+4.08, 22 min wall) showed
# v3 architecture works. v2 reference: 28477675.
#
# Pair with stage_a_amp_pub.sh (B.4-full-amp).

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_standard.sbatch \
  --name stage_a_sup_v3_pub --ntasks 76 --time 12:00:00 \
  --cmd 'tidal sample examples/data/dark_photon_plasma.json \
    --prior "mA2=log_uniform:1e-3:1e3" \
    --prior "deltam=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 256 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1200 \
    --num-repeats 10 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/stage_a_sup_v3_pub'
