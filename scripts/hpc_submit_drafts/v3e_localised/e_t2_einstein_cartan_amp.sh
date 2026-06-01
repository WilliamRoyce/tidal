#!/usr/bin/env bash
# Phase E.T2 amp — Einstein-Cartan minimal (algebraic torsion masses only)
#
# Free couplings (3): alpha1, alpha2, alpha3 — three irreducible torsion mass
# sectors (tentor, trator, axitor). Phase E geometry: localized wavepacket +
# dual-Gaussian B-field, modal-solver compatible.
#
# v3 architecture: tachyon-permissive sampling + soft floor (no Hwang-Noh gate).
# Phase E delta from v3: --ic gaussian (wavepacket) instead of plane-wave;
# dual-Gaussian Bpeak/sigB/zc1/zc2 params instead of uniform B0.
#
# PolyChord budget (ndim=3 landscape):
#   nlive = 25*3 = 75
#   num_repeats = 2*3 = 6
#   max_ndead omitted; precision_criterion 0.01 should converge in <1 hr INTR

set -euo pipefail

# Source the frozen Phase E geometry parameters.
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"

# Pre-flight Boccaletti safety check.
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_t2_amp --ntasks 76 \
  --cmd "tidal sample examples/data/torsion_gertsenshtein_einstein_cartan_e_dual_gaussian.json \
    --prior 'alpha1=arctan_uniform:-89:89' \
    --prior 'alpha2=arctan_uniform:-89:89' \
    --prior 'alpha3=arctan_uniform:-89:89' \
    --likelihood 'P_max:maximize' \
    --baseline-formula 'sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2' \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
    --grid-shape ${N} --bounds=0:${L} --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
    --t-end ${T_END} --snapshots ${SNAPSHOTS} \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 75 --num-repeats 6 \
    --precision-criterion 0.01 \
    --output \${RESULTS_DIR}/e_t2_amp"
