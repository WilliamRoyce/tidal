#!/usr/bin/env bash
# T6-minimal + T6 parity-odd landscape pass.
#
# Phase 1a: T6-minimal amp (76r, 4 params, nlive=100, ~3 min)
# Phase 1b: T6-minimal sup (76r, 4 params, nlive=100, ~3 min)
# Phase 2:  T6 amp (38r) + T6 sup (38r) in parallel (20 params, nlive=500, ~100 min/chain)
#           → Expect SIGTERM at ~60%. Resume with run_landscape_batch_t6_resume.sh.
#
# Run INSIDE an interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   export SLURM_JOB_ID=<jobid> SLURM_NTASKS=76
#   nohup bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_landscape_batch_t6.sh \
#     > /rds/user/wr286/hpc-work/tidal/batch_landscape_t6.log 2>&1 &
#
# After the session: bash scripts/hpc_shuttle.sh pull <jobid>
# If T6 SIGTERM'd: export PREV_JOB_ID=<jobid>; bash run_landscape_batch_t6_resume.sh

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

# Shared physical setup (same grid/IC as D2.x)
PHYS_COMMON="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 10 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# ── Phase 1a: T6-minimal amp landscape (76r, 4 params, nlive=100, ~3 min) ──────
echo "=== [$(date +%H:%M:%S)] Phase 1a: T6-minimal amp landscape (76 ranks) ==="
# shellcheck disable=SC2086
mpirun -n 76 tidal sample examples/data/torsion_gertsenshtein_parity_odd_minimal.json \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "d21=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    ${PHYS_COMMON} \
    --nlive 100 --num-repeats 8 \
    --output "${RESULTS_DIR}/t6_minimal_amp_v3"
echo "=== [$(date +%H:%M:%S)] Phase 1a done ==="
echo ""

# ── Phase 1b: T6-minimal sup landscape (76r, 4 params, nlive=100, ~3 min) ──────
echo "=== [$(date +%H:%M:%S)] Phase 1b: T6-minimal sup landscape (76 ranks) ==="
# shellcheck disable=SC2086
mpirun -n 76 tidal sample examples/data/torsion_gertsenshtein_parity_odd_minimal.json \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "d21=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    ${PHYS_COMMON} \
    --nlive 100 --num-repeats 8 \
    --output "${RESULTS_DIR}/t6_minimal_sup_v3"
echo "=== [$(date +%H:%M:%S)] Phase 1b done ==="
echo ""

# ── Phase 2: T6 amp (38r) + T6 sup (38r) in parallel ─────────────────────────
# ndim=20, nlive=500, num_repeats=40 — est ~100 min each at 38r
# Expect SIGTERM at ~60%; resume with run_landscape_batch_t6_resume.sh
echo "=== [$(date +%H:%M:%S)] Phase 2: T6 amp + T6 sup landscape (38+38 ranks) ==="
echo "    ndim=20, nlive=500, num_repeats=40 — est ~100 min each"
echo "    Expect SIGTERM at ~60%%; resume: export PREV_JOB_ID=${SLURM_JOB_ID}"

T6_DATA="examples/data/torsion_gertsenshtein_parity_odd.json"
T6_PRIORS="--prior beta1=arctan_uniform:-89:89
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

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${T6_DATA} ${PHYS_COMMON} ${T6_PRIORS} \
    --likelihood "P_max:maximize" \
    --nlive 500 --num-repeats 40 \
    --output "${RESULTS_DIR}/t6_amp_v3" \
    > "${RESULTS_DIR}/t6_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${T6_DATA} ${PHYS_COMMON} ${T6_PRIORS} \
    --likelihood "P_max:minimize" \
    --nlive 500 --num-repeats 40 \
    --output "${RESULTS_DIR}/t6_sup_v3" \
    > "${RESULTS_DIR}/t6_sup.log" 2>&1 &
PID_B=$!

echo "    T6 amp pid=${PID_A}  T6 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Done — T6 amp rc=${RC_A}  T6 sup rc=${RC_B} ==="
echo ""
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
echo "=== If T6 SIGTERM'd: export PREV_JOB_ID=${SLURM_JOB_ID}"
echo "===                  bash run_landscape_batch_t6_resume.sh ==="
