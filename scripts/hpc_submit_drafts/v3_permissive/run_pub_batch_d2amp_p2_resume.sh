#!/usr/bin/env bash
# Resume amp pub chains for D2.2, D2.3 (from previous INTR) + D2.1 resume
# (from job 29507332 where it was 63% complete).
#
# Run AFTER run_pub_batch_d2amp_p2.sh SIGTERMs.
#
# Run INSIDE a new interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   export PREV_JOB_ID=<previous INTR jobid from p2.sh>
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_pub_batch_d2amp_p2_resume.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_amp_p2_resume.log 2>&1 &
#
# Phase 1 (~20 min): D2.2 + D2.3 amp pub resume at 38+38 ranks (parallel)
# Phase 2 (~12 min): D2.1 amp pub resume at 76 ranks (from job 29507332)
# Total: ~32 min — leaves plenty of time in a 1hr INTR.
#
# After the session: bash scripts/hpc_shuttle.sh pull <jobid>
# Also pull D2.1 result: bash scripts/hpc_shuttle.sh pull 29507332 --src .../hpc_results/29507332

set -uo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

: "${PREV_JOB_ID:?PREV_JOB_ID must be set to the job ID of run_pub_batch_d2amp_p2.sh}"

echo "=== [$(date +%H:%M:%S)] Setting up PolyChord environment ==="
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_SP}"
tar xf $HOME/venv_site.tar -C "${LOCAL_SP}"

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
# D2.1 amp pub was started in job 29507332 — resume from there
D21_DIR="${TIDAL_ROOT}/hpc_results/29507332/d21_barker_amp_v3_pub"
echo "Resuming D2.2/D2.3 from: ${PREV_RESULTS_DIR}"
echo "Resuming D2.1 from:       ${D21_DIR}"
echo ""

DATA="examples/data/torsion_gertsenshtein_general_nonminimal.json"
PHYS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 10 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# ── Phase 1: Resume D2.2 + D2.3 amp pub (38+38 ranks, parallel) ──────────────
echo "=== [$(date +%H:%M:%S)] Phase 1/2: Resume D2.2 + D2.3 amp pub (38+38 ranks) ==="

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} \
    --param chi=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 400 --num-repeats 16 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/d22_shapiro_amp_v3_pub" \
    > "${PREV_RESULTS_DIR}/d22_amp_pub_resume.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 450 --num-repeats 18 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/d23_full_amp_v3_pub" \
    > "${PREV_RESULTS_DIR}/d23_amp_pub_resume.log" 2>&1 &
PID_B=$!

echo "    D2.2 pid=${PID_A}  D2.3 pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Phase 1 done — D2.2 rc=${RC_A}  D2.3 rc=${RC_B} ==="
echo ""

# ── Phase 2: Resume D2.1 amp pub (76 ranks, from job 29507332) ───────────────
echo "=== [$(date +%H:%M:%S)] Phase 2/2: Resume D2.1 amp pub (76 ranks) ==="

# shellcheck disable=SC2086
mpirun -n 76 tidal sample ${DATA} ${PHYS} \
    --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 300 --num-repeats 12 \
    --read-resume \
    --output "${D21_DIR}"
echo "=== [$(date +%H:%M:%S)] Phase 2 done ==="
echo ""

echo "=== ALL DONE ==="
echo "=== Pull D2.2/D2.3: bash scripts/hpc_shuttle.sh pull ${PREV_JOB_ID} --src ${PREV_RESULTS_DIR} ==="
echo "=== Pull D2.1:       bash scripts/hpc_shuttle.sh pull 29507332 --src ${D21_DIR%/*} ==="
