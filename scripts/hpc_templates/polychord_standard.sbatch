#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH --exclude=cpu-r-10,cpu-r-15

# Standard QOS version of polychord_intr.sbatch (no --qos=intr).
# Use when INTR MaxSubmitPU=1 is saturated or wall time exceeds 1 h.
# All other setup (tarball extraction, MPI modules, staleness check)
# is identical to polychord_intr.sbatch.

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"

# Per-job results directory.  Callers writing ``--output ${RESULTS_DIR}/<subdir>``
# get a stable HPC location that ``hpc_shuttle.sh pull <jobid>`` walks.
export RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
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
The tarball at /home/wr286/venv_site.tar is stale.

Refresh it from the local workstation:

    bash scripts/hpc_install_polychord.sh
    bash scripts/hpc_refresh_venv_tar.sh

Then resubmit this job.
EOF
    exit 3
fi

echo "--- Running tidal command (mpirun -n ${SLURM_NTASKS}) ---"
cd "${TIDAL_ROOT}"
MPIRUN_PREFIX="mpirun -n ${SLURM_NTASKS}"
# Auto-prefix mpirun (PolyChord-specific template — see polychord_intr.sbatch
# for the rationale).
${MPIRUN_PREFIX} {{COMMAND}}

echo "=== DONE ==="
date
