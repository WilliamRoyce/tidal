#!/usr/bin/env bash
# Resume T6 parity-odd amp + sup landscape chains (from SIGTERM in run_landscape_batch_t6.sh).
#
# Run INSIDE a new interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   export PREV_JOB_ID=<jobid from run_landscape_batch_t6.sh>
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_landscape_batch_t6_resume.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_landscape_t6_resume.log 2>&1 &
#
# Resumes T6 amp + T6 sup at 38+38 ranks (remaining ~40% ≈ 40 min total).
# After the session: bash scripts/hpc_shuttle.sh pull ${PREV_JOB_ID} --src ...

set -uo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

: "${PREV_JOB_ID:?PREV_JOB_ID must be set to the job ID of run_landscape_batch_t6.sh}"

echo "=== [$(date +%H:%M:%S)] Setting up PolyChord environment ==="
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"
export TIDAL_MODAL_BACKEND=jax

. /etc/profile.d/modules.sh
module load rhel8/default-icl   2>/dev/null || true
module load intel/2019.3.199    2>/dev/null || true
module load intel-oneapi-mpi    2>/dev/null || true

python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py; print('venv OK')"

cd "${TIDAL_ROOT}"
PREV_RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${PREV_JOB_ID}"
echo "Resuming T6 amp+sup from: ${PREV_RESULTS_DIR}"
echo ""

DATA="examples/data/torsion_gertsenshtein_parity_odd.json"
PHYS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 10 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

PRIORS="--prior beta1=arctan_uniform:-89:89
  --prior beta2=arctan_uniform:-89:89
  --prior beta3=arctan_uniform:-89:89
  --prior xi=log_uniform:1e-3:1e3
  --prior delta1=arctan_uniform:-89:89
  --prior chi=arctan_uniform:-89:89
  --prior zeta1=arctan_uniform:-89:89
  --prior zeta2=arctan_uniform:-89:89
  --prior zeta3=arctan_uniform:-89:89
  --prior d14=arctan_uniform:-89:89
  --prior d15=arctan_uniform:-89:89
  --prior d17=arctan_uniform:-89:89
  --prior d19=arctan_uniform:-89:89
  --prior d20=arctan_uniform:-89:89
  --prior d21=arctan_uniform:-89:89
  --prior zt1=arctan_uniform:-89:89
  --prior zt2=arctan_uniform:-89:89
  --prior zt3=arctan_uniform:-89:89
  --prior zt5=arctan_uniform:-89:89
  --prior zt6=arctan_uniform:-89:89"

echo "=== [$(date +%H:%M:%S)] Resume T6 amp + T6 sup landscape (38+38 ranks) ==="

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} ${PRIORS} \
    --likelihood "P_max:maximize" \
    --nlive 500 --num-repeats 40 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/t6_amp_v3" \
    > "${PREV_RESULTS_DIR}/t6_amp_resume.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} ${PRIORS} \
    --likelihood "P_max:minimize" \
    --nlive 500 --num-repeats 40 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/t6_sup_v3" \
    > "${PREV_RESULTS_DIR}/t6_sup_resume.log" 2>&1 &
PID_B=$!

echo "    T6 amp pid=${PID_A}  T6 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Done — T6 amp rc=${RC_A}  T6 sup rc=${RC_B} ==="
echo ""
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${PREV_JOB_ID} --src ${PREV_RESULTS_DIR} ==="
