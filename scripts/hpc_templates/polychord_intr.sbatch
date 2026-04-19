#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks=1
#SBATCH --time={{TIME}}
#SBATCH --output={{REMOTE_ROOT}}/logs/{{JOB_NAME}}_%j.log
#SBATCH --exclude=cpu-r-10,cpu-r-15

# Reusable template for PolyChord nested-sampling jobs on CSD3.
#
# Preconditions (caller must have done these once, not per-job):
#   1. bash scripts/hpc_shuttle.sh push             # tidal source on HPC
#   2. bash scripts/hpc_install_polychord.sh        # pypolychord built & in venv
#   3. bash scripts/hpc_refresh_venv_tar.sh         # /home/wr286/venv_site.tar regenerated
#
# Per-job steps:
#   1. Extract the prebuilt site-packages tarball from /home (NFS) to
#      /tmp on the compute node, sidestepping /rds Lustre hangs that
#      were observed on some compute nodes (see issue #269).
#   2. Sanity-check that pypolychord is importable from the extracted
#      copy — fail loudly if the tarball is stale.
#   3. Run the caller-supplied {{COMMAND}} (tidal sample ...).

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname) ==="
date

# 1. Extract the pre-built venv tarball from /home to local /tmp.
echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

# 2. Staleness check.  The tarball is expected to contain pypolychord
#    (see Preconditions above).  If it doesn't, the rest of the job
#    will either crash or fall through to a disk-heavy path — surface
#    the root cause immediately instead.
echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic" 2>&1; then
    cat >&2 <<EOF
error: pypolychord or anesthetic missing from the extracted venv.
The tarball at /home/wr286/venv_site.tar is stale.

Refresh it from the local workstation:

    bash scripts/hpc_install_polychord.sh     # (if pypolychord itself needs rebuilding)
    bash scripts/hpc_refresh_venv_tar.sh

Then resubmit this job.
EOF
    exit 3
fi

# 3. Run the caller's tidal command.  Defaults worth remembering:
#      --likelihood-backend memory  (default in v0.33+, ~10x speedup)
#      --num-repeats 2*ndim         (Handley's own run_pypolychord.py)
#      --precision-criterion 0.01   (Cobaya's expensive-likelihood rec)
#      --no-clustering              (safe for unimodal posteriors)
echo "--- Running tidal command ---"
cd "${TIDAL_ROOT}"
{{COMMAND}}

echo "=== DONE ==="
date
