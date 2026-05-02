#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3380
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH --exclude=cpu-r-10,cpu-r-15

# Reusable template for PolyChord nested-sampling jobs on CSD3.
#
# Preconditions (caller must have done these once, not per-job):
#   1. bash scripts/hpc_shuttle.sh push             # tidal source on HPC
#   2. bash scripts/hpc_install_polychord.sh        # pypolychord+mpi4py, MPI=1 build
#   3. bash scripts/hpc_refresh_venv_tar.sh         # /home/wr286/venv_site.tar
#
# Per-job steps:
#   1. Extract the prebuilt site-packages tarball from /home (NFS) to
#      /tmp on the compute node, sidestepping /rds Lustre hangs that
#      were observed on some compute nodes (see issue #269).
#   2. Sanity-check that pypolychord + mpi4py + anesthetic are all
#      importable from the extracted copy — fail loudly if the tarball
#      is stale (#292 MPI path needs mpi4py too).
#   3. Run the caller-supplied command (tidal sample ...) via mpirun.
#      pypolychord auto-detects MPI.COMM_WORLD at import time; with N
#      ranks it parallelizes the slice-sampling step up to ~nlive ranks
#      (arXiv:1506.00171 §3).  {{NTASKS}} is the rendered rank count;
#      the caller controls it via the submit wrapper.
#      NOTE on --cpus-per-task=1: _set_single_thread_blas() already
#      pins each worker to 1 BLAS thread (tidal/cli/_sweep.py),
#      so giving each rank >1 CPU just wastes allocation.

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"

# Per-job results directory.  Callers that want a stable on-HPC output
# location for ``--output`` use ``${RESULTS_DIR}/<subdir>``; the shuttle's
# ``pull <jobid>`` walks the recorded --output path and brings the
# whitelisted artefacts back locally.
export RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"

# Campaign dir: stable path shared across INTR resume runs.
# Set by --campaign NAME in hpc_shuttle.sh submit; empty comment when not used.
# Use ${CAMPAIGN_DIR} in --output so PolyChord finds its .resume checkpoint
# between rounds.  All rounds must pass identical --prior/--n-live/--ntasks.
{{CAMPAIGN_DIR_SETUP}}

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
date

# 1. Extract the pre-built venv tarball from /home to local /tmp.
echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

# MPI toolchain must be loaded inside the job too — sbatch inherits a
# stripped environment on CSD3 and the compute node doesn't auto-load
# site modules.  Needed both for libchord.so's MPI symbols and for
# mpirun itself to find the correct PMI.
module load rhel8/default-icl 2>/dev/null || true
module load intel/2019.3.199 2>/dev/null || true
module load intel-oneapi-mpi 2>/dev/null || true

# 2. Staleness check.  Need pypolychord (always), anesthetic (always),
#    and mpi4py (whenever ntasks>1; harmless when serial).
echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py" 2>&1; then
    cat >&2 <<EOF
error: pypolychord, anesthetic, or mpi4py missing from the extracted venv.
The tarball at /home/wr286/venv_site.tar is stale.

Refresh it from the local workstation:

    bash scripts/hpc_install_polychord.sh     # (rebuilds MPI-enabled PolyChord + installs mpi4py)
    bash scripts/hpc_refresh_venv_tar.sh

Then resubmit this job.
EOF
    exit 3
fi

# 3. Run the caller's tidal command under mpirun.  Defaults worth
#    remembering:
#      --likelihood-backend memory  (default in v0.33+, ~10x speedup)
#      --num-repeats 2*ndim         (Handley's own run_pypolychord.py)
#      --precision-criterion 0.01   (Cobaya's expensive-likelihood rec)
#      --no-clustering              (safe for unimodal posteriors)
#    For ntasks=1, mpirun still works as a pass-through.
echo "--- Running tidal command (mpirun -n ${SLURM_NTASKS}) ---"
cd "${TIDAL_ROOT}"
MPIRUN_PREFIX="mpirun -n ${SLURM_NTASKS}"
# Auto-prefix the caller's command with mpirun so 76 SLURM tasks become
# 76 MPI ranks coordinated by pypolychord.  Previously the template
# defined ``MPIRUN_PREFIX`` but left it to the caller to use it; if the
# caller forgot, SLURM spawned ``ntasks`` independent processes and
# pypolychord's ``MPI_Init_thread`` failed (``PMI2_Job_GetId returned 14``).
# This template is PolyChord-specific so the prefix is unconditional.
${MPIRUN_PREFIX} {{COMMAND}}

echo "=== DONE ==="
date
