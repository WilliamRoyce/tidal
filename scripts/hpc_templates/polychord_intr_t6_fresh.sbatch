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

# T6 parity-odd amp + sup landscape chains — FRESH start with fixed soft floor.
#
# Previous chains (29515407 + 29532548) used SOFT_FLOOR_LOGL=-100 which dragged
# logZ to -100 (63% of 20D prior tachyonic via ε·R̃×F couplings).  The
# PolyChord precision criterion (0.01×|-100|=1 nat) was never reachable.
#
# Fix: --soft-floor-logl -15 → logZ ≈ -1 → threshold 0.01 nat, ~1000 dead
# points, ~40 min.  See issue #375.  No --read-resume (contaminated chains).
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_t6_fresh.sbatch \
#     --ntasks 76 --nodes 1 --name t6-fresh \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== T6 parity-odd fresh landscape (soft-floor-logl=-15, floor fix #375) ==="
date

echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf $HOME/venv_site.tar -C "${LOCAL_SP}"

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

mkdir -p "${RESULTS_DIR}"
cd "${TIDAL_ROOT}"

DATA="examples/data/torsion_gertsenshtein_parity_odd.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 10 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 20 priors — d16/d18/zt4 excluded (zero by symmetry in T6 TOML)
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
echo "    ndim=20, nlive=500, num_repeats=40, soft_floor_logl=-15"
echo "    Est. convergence: ~1000 dead pts ≈ 40 min"

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive 500 --num-repeats 40 \
    --output "${RESULTS_DIR}/t6_amp_v4" \
    > "${RESULTS_DIR}/t6_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive 500 --num-repeats 40 \
    --output "${RESULTS_DIR}/t6_sup_v4" \
    > "${RESULTS_DIR}/t6_sup.log" 2>&1 &
PID_B=$!

echo "    T6 amp pid=${PID_A}  T6 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — T6 amp rc=${RC_A}  T6 sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
