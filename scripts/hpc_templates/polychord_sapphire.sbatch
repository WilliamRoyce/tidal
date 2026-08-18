#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p sapphire
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1500
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH --exclude=cpu-r-11

# Sapphire variant of polychord_intr.sbatch (no --qos=intr).
# Standard-queue PolyChord job on the preferred sapphire partition (112 cores/node).
# Use when wall time exceeds the 1 h INTR limit and multiple-node parallelism is wanted.
#
# Caller {{COMMAND}} is responsible for launching its own mpirun calls — this
# template does NOT auto-prefix mpirun (unlike polychord_intr.sbatch), because
# the typical sapphire use case runs MULTIPLE independent chains in parallel
# inside one allocation (one mpirun per chain) rather than one big chain
# spanning all tasks. Caller uses MPIRUN_PREFIX (provided) or wraps each
# chain command in its own mpirun -n N call.

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"

export RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"

# Campaign dir hook (unused by overnight; reserved for future resume runs).
{{CAMPAIGN_DIR_SETUP}}

echo "=== Job ${SLURM_JOB_ID} on $(hostname), nodes=${SLURM_NNODES} ntasks=${SLURM_NTASKS} ==="
date

# 1. Extract the pre-built venv tarball from /home to local /tmp.
echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf $HOME/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

# 2. Sapphire-specific module setup (rhel8/default-sar per CLAUDE.md).
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-sar
module load python/3.11.0-icl

# MPI toolchain for PolyChord MPI ranks
module load intel/2019.3.199          2>/dev/null || true
module load intel-oneapi-mpi          2>/dev/null || true

# 3. Staleness check (same as polychord_intr.sbatch).
echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py" 2>&1; then
    cat >&2 <<EOF
error: pypolychord, anesthetic, or mpi4py missing from the extracted venv.
The tarball at $HOME/venv_site.tar is stale.

Refresh it from the local workstation:

    bash scripts/hpc_install_polychord.sh
    bash scripts/hpc_refresh_venv_tar.sh

Then resubmit this job.
EOF
    exit 3
fi

# 4. Caller's command. NOT auto-prefixed with mpirun (unlike INTR template).
# Caller must launch its own mpirun calls — this is the parallel-chains pattern.
# MPIRUN_PREFIX is exposed for convenience if the caller wants the "single big
# chain spanning all tasks" pattern.
export MPIRUN_PREFIX="mpirun -n ${SLURM_NTASKS}"
echo "--- Running caller command (NOT auto-mpirun'd; caller controls rank allocation) ---"
cd "${TIDAL_ROOT}"
{{COMMAND}}

echo "=== DONE ==="
date
