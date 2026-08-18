#!/usr/bin/env bash
# Publication-pass amp chains for D2.0–D2.3 (50×ndim nlive).
#
# Run INSIDE an interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_pub_batch_d2amp.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_amp_pub.log 2>&1 &
#
# Strategy: 2-way parallel mpirun (38 ranks each) then next pair.
#   Pair 1: D2.0 amp pub (5p, nlive=250) + D2.1 amp pub (6p, nlive=300) — ~25 min
#   Pair 2: D2.2 amp pub (8p, nlive=400) + D2.3 amp pub (9p, nlive=450) — ~30 min
# Total: ~55 min — fits in 1hr INTR.
#
# All amp chains converge via precision_criterion (no SIGTERM risk).
# After the session: bash scripts/hpc_shuttle.sh pull <jobid>

set -uo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

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
export RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"
echo "RESULTS_DIR=${RESULTS_DIR}"
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

# ── Pair 1: D2.0 amp pub (5p, nlive=250) + D2.1 amp pub (6p, nlive=300) ──────
echo "=== [$(date +%H:%M:%S)] Pair 1: D2.0 amp pub + D2.1 amp pub (38+38 ranks) ==="

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} \
    --param chi=0 --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 250 --num-repeats 10 \
    --output "${RESULTS_DIR}/d20_bahamonde_amp_v3_pub" \
    > "${RESULTS_DIR}/d20_amp_pub.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} \
    --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 300 --num-repeats 12 \
    --output "${RESULTS_DIR}/d21_barker_amp_v3_pub" \
    > "${RESULTS_DIR}/d21_amp_pub.log" 2>&1 &
PID_B=$!

echo "    D2.0 pid=${PID_A}  D2.1 pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Pair 1 done — D2.0 rc=${RC_A}  D2.1 rc=${RC_B} ==="
echo ""

# ── Pair 2: D2.2 amp pub (8p, nlive=400) + D2.3 amp pub (9p, nlive=450) ──────
echo "=== [$(date +%H:%M:%S)] Pair 2: D2.2 amp pub + D2.3 amp pub (38+38 ranks) ==="

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
    --output "${RESULTS_DIR}/d22_shapiro_amp_v3_pub" \
    > "${RESULTS_DIR}/d22_amp_pub.log" 2>&1 &
PID_C=$!

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
    --output "${RESULTS_DIR}/d23_full_amp_v3_pub" \
    > "${RESULTS_DIR}/d23_amp_pub.log" 2>&1 &
PID_D=$!

echo "    D2.2 pid=${PID_C}  D2.3 pid=${PID_D} — waiting..."
wait ${PID_C}; RC_C=$?
wait ${PID_D}; RC_D=$?
echo "=== [$(date +%H:%M:%S)] Pair 2 done — D2.2 rc=${RC_C}  D2.3 rc=${RC_D} ==="
echo ""

echo "=== ALL DONE. Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} --src ${RESULTS_DIR} ==="
