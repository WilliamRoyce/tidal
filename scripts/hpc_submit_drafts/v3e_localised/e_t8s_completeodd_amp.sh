#!/usr/bin/env bash
# Phase E.T8s amp — Complete-odd PGT (already R̃²-free; verified)
#
# IMPORTANT scope: complete-odd Lagrangian has 58 free couplings
# (β,ξ,δ,χ,ζ + d14-21 + zt1-6 + xt1-36). nlive=25*58=1450, num_repeats=2*58=116
# blows past the 1-h INTR limit. This script samples ONLY the SHARED-WITH-T7s
# 18-param subset (β,ξ,δ,ζ,χ-equivalent); the parity-odd-only couplings
# (d14-21, zt1-6, xt1-36) are PINNED to 0 via --param. The full 58D version
# would need a standard-queue job (out of scope for Wave 3 INTR).

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_t8s_amp --ntasks 76 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_complete_odd_e_dual_gaussian.json \
    --prior 'beta1=arctan_uniform:-89:89' \
    --prior 'beta2=arctan_uniform:-89:89' \
    --prior 'beta3=arctan_uniform:-89:89' \
    --prior 'xi=log_uniform:1e-3:1e3' \
    --prior 'delta1=arctan_uniform:-89:89' \
    --prior 'zeta1=arctan_uniform:-89:89' \
    --prior 'zeta2=arctan_uniform:-89:89' \
    --prior 'zeta3=arctan_uniform:-89:89' \
    --prior 'chi1=arctan_uniform:-89:89' \
    --prior 'chi2=arctan_uniform:-89:89' \
    --prior 'chi3=arctan_uniform:-89:89' \
    --prior 'chi4=arctan_uniform:-89:89' \
    --prior 'chi5=arctan_uniform:-89:89' \
    --prior 'chi6=arctan_uniform:-89:89' \
    --prior 'chi7=arctan_uniform:-89:89' \
    --prior 'chi8=arctan_uniform:-89:89' \
    --prior 'chi9=arctan_uniform:-89:89' \
    --prior 'chi10=arctan_uniform:-89:89' \
    --likelihood 'P_max:maximize' \
    --baseline-formula 'sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2' \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
    --grid-shape ${N} --bounds=0:${L} --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
    --t-end ${T_END} --snapshots ${SNAPSHOTS} \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 450 --num-repeats 36 \
    --precision-criterion 0.01 \
    --output \${RESULTS_DIR}/e_t8s_amp"
