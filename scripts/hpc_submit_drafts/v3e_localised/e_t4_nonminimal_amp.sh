#!/usr/bin/env bash
# Phase E.T4 amp — Ricci-EM nonminimal (delta1 R̃_{[μν]} F^{μν})
# Free couplings (4): alpha1, alpha2, alpha3 (algebraic torsion masses) + delta1.
# v3 reference: hpc_results/28520217 (T4 amp v3, log Z=-2.26).

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_t4_amp --ntasks 76 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json \
    --prior 'alpha1=arctan_uniform:-89:89' \
    --prior 'alpha2=arctan_uniform:-89:89' \
    --prior 'alpha3=arctan_uniform:-89:89' \
    --prior 'delta1=arctan_uniform:-89:89' \
    --likelihood 'P_max:maximize' \
    --baseline-formula 'sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2' \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
    --grid-shape ${N} --bounds=0:${L} --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
    --t-end ${T_END} --snapshots ${SNAPSHOTS} \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output \${RESULTS_DIR}/e_t4_amp"
