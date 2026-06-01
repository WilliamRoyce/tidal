#!/usr/bin/env bash
# Phase E.T1 amp — Dark-photon-plasma (PGT torsion trace + Maxwell kinetic mixing + plasma mass)
#
# Free couplings (4): alpha3 (torsion-trace Proca mass), xi (vector-torsion kinetic),
# deltam (F·Ftorsion kinetic mixing), mA2 (plasma photon mass).
#
# Mass relationship: physical m_torsion^2 = +2*alpha3 (stable Proca for alpha3 > 0).
# Per feedback_fv_cdt_equivalence_verified.md, sign convention is post-2026-04-24
# (alpha3 > 0 ↔ FV mT2 > 0, both stable).
#
# PolyChord budget (ndim=4):
#   nlive = 25*4 = 100, num_repeats = 2*4 = 8
#
# v3 reference: hpc_results/28474676 (T1 v5 amp, log Z=-0.073, D_KL=0.024 — null)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_geometry.env"
python "$(git rev-parse --show-toplevel)/scripts/v3e_boccaletti_preflight.py" || exit 1

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name e_t1_amp --ntasks 76 \
  --cmd "tidal sample examples/data/dark_photon_plasma_e_dual_gaussian.json \
    --prior 'alpha3=log_uniform:1e-3:1e3' \
    --prior 'xi=log_uniform:1e-3:1e3' \
    --prior 'deltam=arctan_uniform:-89:89' \
    --prior 'mA2=log_uniform:1e-3:1e3' \
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
    --output \${RESULTS_DIR}/e_t1_amp"
