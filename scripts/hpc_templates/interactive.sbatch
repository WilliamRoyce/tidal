#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p sapphire
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --time={{TIME}}
#SBATCH --exclude=cpu-r-11
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out

set -euo pipefail

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-sar
module load python/3.11.0-icl

echo "host=$(hostname) cpus=${SLURM_CPUS_ON_NODE:-?} ntasks=${SLURM_NTASKS:-?} job=${SLURM_JOB_ID}"
echo "started=$(date -Iseconds)"

cd {{REMOTE_ROOT}}
source .venv/bin/activate

export RESULTS_DIR="{{REMOTE_ROOT}}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"

echo "RESULTS_DIR=${RESULTS_DIR}"
echo "=== INTERACTIVE SESSION READY ==="
echo "=== node=$(hostname) ==="
echo "=== attach: bash scripts/hpc_shuttle.sh attach ${SLURM_JOB_ID} ==="

# Print clean finish message even when SLURM kills the sleep via SIGTERM at walltime.
trap 'echo "finished=$(date -Iseconds)"' EXIT

# Keep allocation alive until SLURM kills it at walltime, or scancel.
# sleep infinity is supported on RHEL8 (GNU coreutils); SLURM sends SIGTERM at walltime.
sleep infinity
