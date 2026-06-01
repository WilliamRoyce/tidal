#!/bin/bash
#SBATCH -J atlas-{{THEORY_TAG}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1500
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out

# Phase E atlas — standard-queue sapphire single-node parallel chains for one theory.
#
# Reuses the validated run_atlas_slot.sh runner from the INTR slot pattern.
# Differences vs INTR:
#   - Standard queue (no --qos=intr): longer walltime, queue wait possible
#   - Single theory per sbatch (no packing): simpler scheduling, cleaner failure mode
#   - TIDAL_FRESH=1: no resume contention, clean start
#   - Sapphire partition (preferred CSD3 CPU class per CLAUDE.md)
#
# Substitution variables (rendered by submit wrapper):
#   THEORY_TAG     — t5_1 / t5_2 / t5_3
#   ACCOUNT        — billing account (DiRAC > SL2 > SL3 via resolve-account)
#   NTASKS         — total cores requested (≤112 for single sapphire node)
#   TIME           — walltime budget (HH:MM:SS)
#   REMOTE_ROOT    — /rds/user/wr286/hpc-work/tidal
#
# Per-theory ranks/face is encoded in run_atlas_slot.sh's case statement.
# Empirical sizing (from today's T5.0 runs at 11 ranks/face → ~80 min total):
#   T5.1: 12 chains × 9 ranks = 108 cores, ~3-4h wall  (--time 06:00:00 cushion)
#   T5.2: 16 chains × 7 ranks = 112 cores, ~5-6h wall  (--time 09:00:00)
#   T5.3: 18 chains × 6 ranks = 108 cores, ~8-10h wall (--time 12:00:00)
#
# OOM defense: uniform-field grid is tiny (~100 KB/rank state). 1500 MB/cpu
# is overkill but matches the validated INTR runs.

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
THEORY_TAG="{{THEORY_TAG}}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== Atlas standard-queue: theory=${THEORY_TAG} ==="
date

# 1. Extract pre-built site-packages tarball (avoids Lustre import hangs).
echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

# 2. MPI + module setup (sbatch inherits stripped env).
. /etc/profile.d/modules.sh 2>/dev/null || true
module load rhel8/default-sar 2>/dev/null || true
module load python/3.11.0-icl 2>/dev/null || true
module load intel/2019.3.199 2>/dev/null || true
module load intel-oneapi-mpi 2>/dev/null || true

# 3. Staleness check on extracted venv.
echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py" 2>&1; then
    cat >&2 <<EOF
error: pypolychord, anesthetic, or mpi4py missing from extracted venv.
Refresh:  bash scripts/hpc_install_polychord.sh && bash scripts/hpc_refresh_venv_tar.sh
EOF
    exit 3
fi

mkdir -p "${RESULTS_DIR}"
cd "${TIDAL_ROOT}"

# 4. Export RESULTS_DIR so run_atlas_slot.sh picks it up.
#    TIDAL_FRESH=1: defensive clean start (no stale .resume files).
#    TIDAL_SEQUENTIAL unset (default parallel): all faces launch concurrently.
export RESULTS_DIR
export SLURM_JOB_ID
export TIDAL_FRESH=1

echo "--- Launching run_atlas_slot.sh ${THEORY_TAG} (parallel, fresh) ---"
echo "    Per-theory ndim/faces/ranks from run_atlas_slot.sh case statement."
echo "    Walltime budget: {{TIME}}"
echo

bash scripts/hpc_submit_drafts/v3_atlas/run_atlas_slot.sh "${THEORY_TAG}"
RC=$?

echo
echo "=== [$(date +%H:%M:%S)] DONE — theory=${THEORY_TAG} rc=${RC} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
exit ${RC}
