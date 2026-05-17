#!/usr/bin/env bash
# Second interactive batch: D2.3 sup + parallel mpirun experiment.
#
# Phase 1: D2.3 sup at full 76 ranks (the one theory still missing).
# Phase 2: Test parallel mpirun on one sapphire node — run two publication-
#          resolution reruns simultaneously with 38 MPI ranks each.
#          If this works, future campaigns can pair theories for ~2× throughput.
#
# Run via:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   <ssh to compute node, set SLURM_JOB_ID, SLURM_NTASKS>
#   bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_interactive_batch_2.sh

set -uo pipefail   # NOT -e: keep going if Phase 2 mpirun parallel fails

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

# ── Phase 1: D2.3 Full sup (9p) at full 76 ranks ────────────────────────────
echo "=== [$(date +%H:%M:%S)] Phase 1/2: D2.3 Full sup (76 ranks, nlive=225) ==="
# shellcheck disable=SC2086
mpirun -n 76 tidal sample ${DATA} ${PHYS} \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --nlive 225 --num-repeats 18 \
    --output "${RESULTS_DIR}/d23_full_sup_v3"
echo "=== [$(date +%H:%M:%S)] Phase 1 done ==="
echo ""

# ── Phase 2: Parallel mpirun test — two publication-quality reruns ──────────
#
# Each mpirun gets 38 ranks. If IntelMPI handles dual-mpirun on one node
# cleanly, both finish in ~2x solo time. Otherwise we see PMI errors here
# and learn that sequential is the only safe pattern.
echo "=== [$(date +%H:%M:%S)] Phase 2/2: Parallel mpirun test (38+38 ranks) ==="
echo "    A: D2.0 sup publication (nlive=250, num_repeats=10)"
echo "    B: D2.1 sup publication (nlive=300, num_repeats=12)"

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS} \
    --param chi=0 --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --nlive 250 --num-repeats 10 \
    --output "${RESULTS_DIR}/d20_bahamonde_sup_v3_pub" \
    > "${RESULTS_DIR}/d20_bahamonde_sup_v3_pub.log" 2>&1 &
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
    --likelihood "P_max:minimize" \
    --nlive 300 --num-repeats 12 \
    --output "${RESULTS_DIR}/d21_barker_sup_v3_pub" \
    > "${RESULTS_DIR}/d21_barker_sup_v3_pub.log" 2>&1 &
PID_B=$!

echo "    A pid=${PID_A}  B pid=${PID_B} — waiting on both..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?
echo "=== [$(date +%H:%M:%S)] Phase 2 done — A rc=${RC_A}  B rc=${RC_B} ==="
echo ""
echo "=== ALL DONE. Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} --src ${RESULTS_DIR} ==="
