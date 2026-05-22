#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1500
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH --exclude=cpu-r-10,cpu-r-15

# Resume T6 parity-odd amp + sup landscape chains (38+38 ranks in parallel).
#
# Picks up PolyChord .resume checkpoints from job 29515407 (~60% complete,
# nlive=500, ndim=20).  Estimated remaining wall: ~40 min.
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_t6_resume.sbatch \
#     --ntasks 76 --nodes 1 --name t6-resume \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull 29515407

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
PREV_JOB_ID="29515407"
PREV_RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${PREV_JOB_ID}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== Resuming T6-full amp+sup from ${PREV_RESULTS_DIR} ==="
date

echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

module load rhel8/default-icl 2>/dev/null || true
module load intel/2019.3.199 2>/dev/null || true
module load intel-oneapi-mpi 2>/dev/null || true

echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py" 2>&1; then
    cat >&2 <<EOF
error: pypolychord, anesthetic, or mpi4py missing from the extracted venv.
Refresh: bash scripts/hpc_install_polychord.sh && bash scripts/hpc_refresh_venv_tar.sh
EOF
    exit 3
fi

cd "${TIDAL_ROOT}"

DATA="examples/data/torsion_gertsenshtein_parity_odd.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 10 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

PRIOR_ARGS="--prior beta1=arctan_uniform:-89:89
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

echo "--- Launching T6 amp (38r) + T6 sup (38r) in parallel ---"

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive 500 --num-repeats 40 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/t6_amp_v3" \
    > "${PREV_RESULTS_DIR}/t6_amp_resume2.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive 500 --num-repeats 40 \
    --read-resume \
    --output "${PREV_RESULTS_DIR}/t6_sup_v3" \
    > "${PREV_RESULTS_DIR}/t6_sup_resume2.log" 2>&1 &
PID_B=$!

echo "    T6 amp pid=${PID_A}  T6 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — T6 amp rc=${RC_A}  T6 sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${PREV_JOB_ID} ==="
