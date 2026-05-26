#!/usr/bin/env bash
# Phase E overnight LANDSCAPE-quality run for E.T7s ONLY (complete-even, 18D).
# Separate job from the Wave-1+2 overnight because T7s' 18-dim prior
# (nlive=450, num_repeats=36) is the largest budget by far, and combining
# it with the smaller theories pushes total cores above the SL3 QOS
# cpu=448 limit when using 32 ranks per chain.
#
# At 128 ranks per chain × 2 chains = 256 cores total → fits cpu=448 limit
# easily; estimated wall time ~5-6h at LANDSCAPE quality.

set -euo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"
GEOMETRY_ENV="${TIDAL_ROOT}/scripts/hpc_submit_drafts/v3e_localised/_geometry.env"

# shellcheck disable=SC1090
source "${GEOMETRY_ENV}"

cd "${TIDAL_ROOT}"
: "${RESULTS_DIR:?RESULTS_DIR must be set by sbatch template}"

COMMON_PHYS="\
  --baseline-formula sin(kappa*Bpeak*sigB*sqrt(2*pi)/2)**2 \
  --soft-floor-noise 1.0 \
  --param kappa=1.0 --param Bpeak=${BPEAK} --param sigB=${SIGB} --param zc1=${ZC1} --param zc2=${ZC2} \
  --grid-shape ${N} --bounds=0:${L} --periodic \
  --ic gaussian --ic-component h_5 --ic-amplitude ${H0} --ic-width ${SIGMA_W} --ic-center ${X_C} --ic-wavevector ${K_CARRIER} \
  --t-end ${T_END} --snapshots ${SNAPSHOTS} \
  --measure conversion,peak_conversion --source h_5 --target a_1 \
  --method nested --sampler polychord --precision-criterion 0.01"

# E.T7s priors: 18 BSM couplings (Bahamonde 5 + zeta1-3 + chi1-10)
PRIORS_T7S="\
  --prior beta1=arctan_uniform:-89:89 --prior beta2=arctan_uniform:-89:89 \
  --prior beta3=arctan_uniform:-89:89 --prior xi=log_uniform:1e-3:1e3 \
  --prior delta1=arctan_uniform:-89:89 \
  --prior zeta1=arctan_uniform:-89:89 --prior zeta2=arctan_uniform:-89:89 \
  --prior zeta3=arctan_uniform:-89:89 \
  --prior chi1=arctan_uniform:-89:89 --prior chi2=arctan_uniform:-89:89 \
  --prior chi3=arctan_uniform:-89:89 --prior chi4=arctan_uniform:-89:89 \
  --prior chi5=arctan_uniform:-89:89 --prior chi6=arctan_uniform:-89:89 \
  --prior chi7=arctan_uniform:-89:89 --prior chi8=arctan_uniform:-89:89 \
  --prior chi9=arctan_uniform:-89:89 --prior chi10=arctan_uniform:-89:89"

JSON_T7S="examples/data/torsion_gertsenshtein_complete_even_e_dual_gaussian.json"

echo "=== [$(date +%H:%M:%S)] E.T7s LANDSCAPE: 2 parallel chains (amp+sup) ==="

# shellcheck disable=SC2086
srun --exclusive --ntasks=128 --cpus-per-task=1 --mem-per-cpu=1500 --mpi=pmi2 \
  tidal sample "${JSON_T7S}" ${PRIORS_T7S} \
  --likelihood "P_max:maximize" \
  ${COMMON_PHYS} \
  --nlive 450 --num-repeats 36 \
  --output "${RESULTS_DIR}/e_t7s_amp" \
  > "${RESULTS_DIR}/e_t7s_amp.log" 2>&1 &
PID_AMP=$!

# shellcheck disable=SC2086
srun --exclusive --ntasks=128 --cpus-per-task=1 --mem-per-cpu=1500 --mpi=pmi2 \
  tidal sample "${JSON_T7S}" ${PRIORS_T7S} \
  --likelihood "P_max:minimize" \
  ${COMMON_PHYS} \
  --nlive 450 --num-repeats 36 \
  --output "${RESULTS_DIR}/e_t7s_sup" \
  > "${RESULTS_DIR}/e_t7s_sup.log" 2>&1 &
PID_SUP=$!

echo "=== [$(date +%H:%M:%S)] Both T7s chains launched: amp=${PID_AMP} sup=${PID_SUP} ==="

amp_status=0; sup_status=0
wait ${PID_AMP} || amp_status=$?
wait ${PID_SUP} || sup_status=$?

echo "=== [$(date +%H:%M:%S)] T7s DONE — amp_status=${amp_status} sup_status=${sup_status} ==="
date
