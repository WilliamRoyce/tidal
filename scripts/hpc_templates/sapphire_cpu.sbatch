#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p sapphire
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

{{COMMAND}}

echo "finished=$(date -Iseconds)"
