#!/usr/bin/env bash
# Resume D1 amp + D1 sup landscape chains (from previous INTR SIGTERM).
#
# Run AFTER run_landscape_batch_stagea_d1.sh SIGTERMs on the D1 parallel pair.
#
# Run INSIDE a new interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   export PREV_JOB_ID=<jobid from run_landscape_batch_stagea_d1.sh>
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_landscape_batch_d1_resume.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_landscape_d1_resume.log 2>&1 &
#
# Resumes D1 amp + D1 sup at 38+38 ranks (remaining ~10% each ≈ 6 min total).
# After the session: bash scripts/hpc_shuttle.sh pull ${PREV_JOB_ID} --src ...

set -uo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

: "${PREV_JOB_ID:?PREV_JOB_ID must be set to the job ID of run_landscape_batch_stagea_d1.sh}"

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
echo "Resuming D1 amp+sup from: ${PREV_RESULTS_DIR}"
echo ""

echo "=== [$(date +%H:%M:%S)] Resume D1 amp + D1 sup landscape (38+38 ranks) ==="

mpirun -n 38 tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/d1_amp_v3" \
    > "${PREV_RESULTS_DIR}/d1_amp_resume.log" 2>&1 &
PID_A=$!

mpirun -n 38 tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/d1_sup_v3" \
    > "${PREV_RESULTS_DIR}/d1_sup_resume.log" 2>&1 &
PID_B=$!

echo "    D1 amp pid=${PID_A}  D1 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Done — D1 amp rc=${RC_A}  D1 sup rc=${RC_B} ==="
echo ""
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${PREV_JOB_ID} --src ${PREV_RESULTS_DIR} ==="
