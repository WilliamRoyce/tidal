#!/usr/bin/env bash
# Phase E.NP amp — Non-propagating-torsion variant
# Reuses E.T5.3 JSON (general_nonminimal) with xi pinned to 0 → vector torsion
# becomes purely algebraic (no kinetic; constraint-solved). All other couplings
# free: beta1-3 (algebraic torsion masses) + delta1 + chi + zeta1-3. ndim=8.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_np_amp --ntasks 76 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json \
    --prior 'beta1=arctan_uniform:-89:89' \
    --prior 'beta2=arctan_uniform:-89:89' \
    --prior 'beta3=arctan_uniform:-89:89' \
    --prior 'delta1=arctan_uniform:-89:89' \
    --prior 'chi=arctan_uniform:-89:89' \
    --prior 'zeta1=arctan_uniform:-89:89' \
    --prior 'zeta2=arctan_uniform:-89:89' \
    --prior 'zeta3=arctan_uniform:-89:89' \
    --likelihood 'P_max:maximize' \
    --baseline-formula 'sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2' \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} --param xi=0.0 \
    --grid-shape ${N} --bounds=0:${L} --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
    --t-end ${T_END} --snapshots ${SNAPSHOTS} \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 200 --num-repeats 16 \
    --precision-criterion 0.01 \
    --output \${RESULTS_DIR}/e_np_amp"
