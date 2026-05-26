#!/usr/bin/env bash
# Smoke test for the srun + pypolychord launch pattern before committing
# the overnight job. Validates that parallel srun calls within one SLURM
# allocation properly carve resources and that pypolychord initializes
# MPI correctly under srun --mpi=pmi2.
#
# Tiny chain: nlive=8, num_repeats=2, precision=0.1 — converges in
# seconds (~10 likelihood evals per ndead × tiny ndead count).
# Single theory (E.T4) × amp + sup = 2 parallel chains × 16 ranks each
# = 32 cores out of one INTR sapphire allocation (112 cores).

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
  --method nested --sampler polychord --precision-criterion 0.1"

PRIORS_T4="--prior alpha1=arctan_uniform:-89:89 --prior alpha2=arctan_uniform:-89:89 \
  --prior alpha3=arctan_uniform:-89:89 --prior delta1=arctan_uniform:-89:89"

JSON_T4="examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json"

echo "=== [$(date +%H:%M:%S)] SMOKE: 2 parallel chains (T4 amp + sup) via srun ==="

# shellcheck disable=SC2086
srun --exclusive --ntasks=16 --cpus-per-task=1 --mem-per-cpu=1500 --mpi=pmi2 \
  tidal sample "${JSON_T4}" ${PRIORS_T4} --likelihood "P_max:maximize" \
  ${COMMON_PHYS} --nlive 8 --num-repeats 2 \
  --output "${RESULTS_DIR}/smoke_t4_amp" > "${RESULTS_DIR}/smoke_t4_amp.log" 2>&1 &
PID_AMP=$!

# shellcheck disable=SC2086
srun --exclusive --ntasks=16 --cpus-per-task=1 --mem-per-cpu=1500 --mpi=pmi2 \
  tidal sample "${JSON_T4}" ${PRIORS_T4} --likelihood "P_max:minimize" \
  ${COMMON_PHYS} --nlive 8 --num-repeats 2 \
  --output "${RESULTS_DIR}/smoke_t4_sup" > "${RESULTS_DIR}/smoke_t4_sup.log" 2>&1 &
PID_SUP=$!

echo "=== [$(date +%H:%M:%S)] Both chains launched (amp=${PID_AMP} sup=${PID_SUP}) — waiting ==="

amp_status=0; sup_status=0
wait ${PID_AMP} || amp_status=$?
wait ${PID_SUP} || sup_status=$?

echo
echo "=== [$(date +%H:%M:%S)] SMOKE DONE — amp_status=${amp_status} sup_status=${sup_status} ==="
echo
echo "=== amp log tail ==="
tail -15 "${RESULTS_DIR}/smoke_t4_amp.log" 2>&1
echo
echo "=== sup log tail ==="
tail -15 "${RESULTS_DIR}/smoke_t4_sup.log" 2>&1
echo
echo "=== artefacts ==="
ls -la "${RESULTS_DIR}/smoke_t4_amp/" "${RESULTS_DIR}/smoke_t4_sup/" 2>&1 | tail -20

if [ "${amp_status}" -eq 0 ] && [ "${sup_status}" -eq 0 ] && \
   [ -f "${RESULTS_DIR}/smoke_t4_amp/inference.json" ] && \
   [ -f "${RESULTS_DIR}/smoke_t4_sup/inference.json" ]; then
  echo
  echo "=== SMOKE PASS — srun + pypolychord parallel pattern WORKS ==="
else
  echo
  echo "=== SMOKE FAIL — see logs above ==="
  exit 1
fi
