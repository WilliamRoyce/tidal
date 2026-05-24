#!/usr/bin/env bash
# Phase E.NP sup — companion to e_np_algebraic_amp.sh (minimise P_max).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_np_sup --ntasks 76 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json \
    --prior 'beta1=arctan_uniform:-89:89' \
    --prior 'beta2=arctan_uniform:-89:89' \
    --prior 'beta3=arctan_uniform:-89:89' \
    --prior 'delta1=arctan_uniform:-89:89' \
    --prior 'chi=arctan_uniform:-89:89' \
    --prior 'zeta1=arctan_uniform:-89:89' \
    --prior 'zeta2=arctan_uniform:-89:89' \
    --prior 'zeta3=arctan_uniform:-89:89' \
    --likelihood 'P_max:minimize' \
    --baseline-formula 'sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2' \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} --param xi=0.0 \
    --grid-shape ${N} --bounds=0:${L} --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
    --t-end ${T_END} --snapshots ${SNAPSHOTS} \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 200 --num-repeats 16 \
    --precision-criterion 0.01 \
    --output \${RESULTS_DIR}/e_np_sup"
