#!/usr/bin/env bash
# Phase E atlas — submit one (or more) theories to the STANDARD-queue sapphire.
#
# Usage:
#   bash submit_atlas_standard.sh <theory_tag> [<theory_tag> ...]
#
# Examples:
#   bash submit_atlas_standard.sh t5_1                 # single theory, 6h walltime
#   bash submit_atlas_standard.sh t5_1 t5_2 t5_3       # three sbatches in queue
#
# Each theory becomes a SEPARATE sbatch (no packing) so:
#   - Scheduling is faster (smaller jobs schedule sooner)
#   - One theory's failure doesn't kill the others
#   - Per-theory pull + render is independent
#
# Per-theory sizing (validated by today's INTR runs + extrapolation):
#   t5_1: 108 cores (12×9 ranks), 6h walltime
#   t5_2: 112 cores (16×7 ranks), 9h walltime
#   t5_3: 108 cores (18×6 ranks), 12h walltime
#
# Pre-flight checks (do NOT skip these — they prevent silent failures):
#   1. resolve-account picks DiRAC > SL2 > SL3 (no silent default to SL3)
#   2. tarball staleness verified (sbatch fails fast if pypolychord missing)
#   3. ~30s after each submit, sacct check catches FAILED/CANCELLED early
#
# After submit:
#   bash scripts/hpc_shuttle.sh wait  <jobid>     # blocks until SLURM marks done
#   bash scripts/hpc_shuttle.sh pull  <jobid>     # pulls inference + chains
#   uv run tidal plot hpc_results/<jobid>/atlas_<theory>/ --type atlas \
#       --output hpc_results/<jobid>/figures/atlas_<theory>.pdf

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <theory_tag> [<theory_tag> ...]"
  echo "  theory tags: t5_1 t5_2 t5_3 (T5.0 already complete; T4/T1 done in S1)"
  exit 1
fi

REPO_ROOT="$(git rev-parse --show-toplevel)"
SHUTTLE="${REPO_ROOT}/scripts/hpc_shuttle.sh"
TEMPLATE="${REPO_ROOT}/scripts/hpc_submit_drafts/v3_atlas/sbatch_standard_atlas.sbatch"

if [[ ! -f "${TEMPLATE}" ]]; then
  echo "ERROR: missing template: ${TEMPLATE}" >&2
  exit 1
fi

# Per-theory resource map (must match run_atlas_slot.sh case statement)
sizing_for() {
  case "$1" in
    t5_1) NTASKS=108; TIME="06:00:00" ;;
    t5_2) NTASKS=112; TIME="09:00:00" ;;
    t5_3) NTASKS=108; TIME="12:00:00" ;;
    t5_0) NTASKS=110; TIME="02:00:00" ;;  # if needed for re-runs
    *)
      echo "ERROR: unsupported theory tag for standard queue: $1" >&2
      echo "       Supported: t5_0 t5_1 t5_2 t5_3" >&2
      exit 2
      ;;
  esac
}

# 1. Resolve account once (DiRAC > SL2 > SL3 per hpc_shuttle.sh)
ACCOUNT_LINE="$(bash "${SHUTTLE}" resolve-account 2>&1 | tail -1)"
ACCOUNT="${ACCOUNT_LINE}"
if [[ -z "${ACCOUNT}" || "${ACCOUNT}" =~ [[:space:]] ]]; then
  echo "ERROR: failed to resolve account cleanly: '${ACCOUNT_LINE}'" >&2
  bash "${SHUTTLE}" resolve-account 2>&1 | tail -5
  exit 1
fi
echo "=== Resolved account: ${ACCOUNT} ==="
echo

REMOTE_ROOT="/rds/user/wr286/hpc-work/tidal"
JOBIDS=()
THEORIES=()

for tag in "$@"; do
  sizing_for "${tag}"
  echo "=== Submitting atlas-${tag}: ntasks=${NTASKS}, time=${TIME} ==="

  # Render the sbatch from template
  RENDERED="/tmp/sbatch_atlas_${tag}_$$.sbatch"
  sed \
    -e "s|{{THEORY_TAG}}|${tag}|g" \
    -e "s|{{ACCOUNT}}|${ACCOUNT}|g" \
    -e "s|{{NTASKS}}|${NTASKS}|g" \
    -e "s|{{TIME}}|${TIME}|g" \
    -e "s|{{REMOTE_ROOT}}|${REMOTE_ROOT}|g" \
    "${TEMPLATE}" > "${RENDERED}"

  # Upload + submit via shuttle (uses validated path)
  REMOTE_SBATCH="${REMOTE_ROOT}/slurm_logs/sbatch_atlas_${tag}.sbatch"
  ssh csd3 "mkdir -p ${REMOTE_ROOT}/slurm_logs"
  scp -q "${RENDERED}" "csd3:${REMOTE_SBATCH}"
  jobid="$(ssh csd3 "cd ${REMOTE_ROOT} && sbatch ${REMOTE_SBATCH}" | awk '{print $NF}')"
  rm -f "${RENDERED}"

  if [[ ! "${jobid}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: sbatch did not return a numeric jobid: '${jobid}'" >&2
    exit 1
  fi
  echo "  Submitted: jobid=${jobid}"
  JOBIDS+=("${jobid}")
  THEORIES+=("${tag}")
  echo
done

# 2. Fast-fail check: sleep 30s, then sacct check ALL jobids
echo "=== Fast-fail check: sleeping 30s then sacct ==="
sleep 30
for i in "${!JOBIDS[@]}"; do
  jobid="${JOBIDS[$i]}"
  tag="${THEORIES[$i]}"
  state="$(ssh csd3 "sacct -j ${jobid} --format=State -P --noheader 2>/dev/null | head -1")"
  case "${state}" in
    FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY)
      echo "  WARNING: ${tag} (jobid ${jobid}) state=${state} — investigate immediately"
      ;;
    PENDING|RUNNING|"")
      echo "  OK: ${tag} (jobid ${jobid}) state=${state:-PENDING}"
      ;;
    *)
      echo "  UNKNOWN: ${tag} (jobid ${jobid}) state=${state}"
      ;;
  esac
done

echo
echo "=== Submitted ${#JOBIDS[@]} sbatches ==="
for i in "${!JOBIDS[@]}"; do
  echo "  ${THEORIES[$i]}: jobid=${JOBIDS[$i]}"
done
echo
echo "Monitor:   ssh csd3 \"squeue -j $(IFS=,; echo "${JOBIDS[*]}")\""
echo "Pull each: bash scripts/hpc_shuttle.sh pull <jobid>"
