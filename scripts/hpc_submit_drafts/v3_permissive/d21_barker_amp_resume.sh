#!/usr/bin/env bash
# v3 D2.1 Barker amp — RESUME from 29256858 (TIMEOUT expected: 27 clusters, slow convergence)
#
# 29256858 found 27 posterior clusters — multi-modal posterior requires more iterations.
# logZ = 9.56 ± 0.049 at 21 min; precision_criterion=0.01 not yet met.
# Resume from checkpoint; output written back to SAME directory (pull uses: hpc_shuttle.sh pull 29256858).

set -euo pipefail

PREV_OUTPUT=/rds/user/wr286/hpc-work/tidal/hpc_results/29256858/d21_barker_amp_v3

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d21_barker_amp_v3r --ntasks 76 --time 01:00:00 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
    --param kappa=1.0 --param B0=0.01 \
    --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior 'beta1=arctan_uniform:-89:89' \
    --prior 'beta2=arctan_uniform:-89:89' \
    --prior 'beta3=arctan_uniform:-89:89' \
    --prior 'xi=log_uniform:1e-3:1e3' \
    --prior 'delta1=arctan_uniform:-89:89' \
    --prior 'chi=arctan_uniform:-89:89' \
    --likelihood 'P_max:maximize' \
    --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
    --soft-floor-noise 1.0 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1800 \
    --num-repeats 5 --precision-criterion 0.01 \
    --output ${PREV_OUTPUT} \
    --read-resume"
