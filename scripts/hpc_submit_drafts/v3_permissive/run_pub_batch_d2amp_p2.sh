#!/usr/bin/env bash
# Publication-pass amp chains for D2.2 and D2.3 (50×ndim nlive) — Pair 2 only.
#
# D2.0 and D2.1 amp pub already captured in job 29507332 (D2.0 complete,
# D2.1 partial — resume via run_pub_batch_d2amp_p2_resume.sh).
#
# Run INSIDE an interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_pub_batch_d2amp_p2.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_amp_p2_pub.log 2>&1 &
#
# Strategy: 2-way parallel mpirun (38 ranks each).
#   D2.2 amp pub (8p, nlive=400): est ~80 min at 38 ranks → SIGTERM at ~75%
#   D2.3 amp pub (9p, nlive=450): est ~80 min at 38 ranks → SIGTERM at ~75%
#
# After SIGTERM: pull the job then resume with run_pub_batch_d2amp_p2_resume.sh
# (sets PREV_JOB_ID=<this jobid> before running).
#
# After the session: bash scripts/hpc_shuttle.sh pull <jobid>

set -uo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

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

echo "=== [$(date +%H:%M:%S)] D2.2 amp pub + D2.3 amp pub (38+38 ranks) ==="

# ── D2.2 Shapiro amp pub (8p: β₁₋₃, ξ, δ₁, ζ₁₋₃) ──────────────────────────
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
PID_A=$!

# ── D2.3 Full amp pub (9p: all couplings) ────────────────────────────────────
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
PID_B=$!

echo "    D2.2 pid=${PID_A}  D2.3 pid=${PID_B} — waiting (expect SIGTERM at 1hr ~75% done)..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Done — D2.2 rc=${RC_A}  D2.3 rc=${RC_B} ==="
echo ""
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
echo "=== Resume: export PREV_JOB_ID=${SLURM_JOB_ID}"
echo "===         bash run_pub_batch_d2amp_p2_resume.sh ==="
