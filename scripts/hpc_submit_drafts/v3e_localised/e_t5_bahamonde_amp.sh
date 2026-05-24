#!/usr/bin/env bash
# Phase E.T5 amp — Bahamonde 5D subspace of general_nonminimal
# Free couplings (5): beta1, beta2, beta3, xi, delta1
# Pin chi=0 and zeta1-3=0 via fixed-prior at 0 → restricts to Bahamonde subspace.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_t5_amp --ntasks 76 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json \
    --prior 'beta1=arctan_uniform:-89:89' \
    --prior 'beta2=arctan_uniform:-89:89' \
    --prior 'beta3=arctan_uniform:-89:89' \
    --prior 'xi=log_uniform:1e-3:1e3' \
    --prior 'delta1=arctan_uniform:-89:89' \
    --likelihood 'P_max:maximize' \
    --baseline-formula 'sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2' \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
    --param chi=0.0 --param zeta1=0.0 --param zeta2=0.0 --param zeta3=0.0 \
    --grid-shape ${N} --bounds=0:${L} --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
    --t-end ${T_END} --snapshots ${SNAPSHOTS} \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 125 --num-repeats 10 \
    --precision-criterion 0.01 \
    --output \${RESULTS_DIR}/e_t5_amp"
