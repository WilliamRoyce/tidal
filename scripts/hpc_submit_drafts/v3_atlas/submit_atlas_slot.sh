#!/usr/bin/env bash
# Phase E Atlas — local wrapper that books an INTR slot, waits for READY,
# then ssh-launches run_atlas_slot.sh on the compute node.
#
# Usage:
#   bash submit_atlas_slot.sh <theory_tag> [<theory_tag> ...]
#
# Examples:
#   Phase α S1:  bash submit_atlas_slot.sh t4 t1
#   Phase α S2:  bash submit_atlas_slot.sh t5_0 t5_1
#
# Returns the INTR slot jobid + node hostname; subsequent monitoring can use
#   bash scripts/hpc_shuttle.sh wait <jobid>
#   bash scripts/hpc_shuttle.sh pull <jobid>
# and BATCH_DONE marker on /rds/.../hpc_results/<jobid>/BATCH_DONE.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <theory_tag> [<theory_tag> ...]"
  exit 1
fi

THEORIES=("$@")
THEORIES_STR="${THEORIES[*]}"
REPO_ROOT="$(git rev-parse --show-toplevel)"
SHUTTLE="${REPO_ROOT}/scripts/hpc_shuttle.sh"

# 1. Book INTR slot
echo "=== Booking INTR sapphire slot for atlas theories: ${THEORIES_STR} ==="
jobid_line="$(bash "${SHUTTLE}" interactive 2>&1 | grep -E '^[0-9]+$' | head -1)"
SLOT_JOBID="${jobid_line}"

if [[ -z "${SLOT_JOBID}" || ! "${SLOT_JOBID}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: failed to parse jobid from interactive command output" >&2
  bash "${SHUTTLE}" interactive 2>&1 | tail -5
  exit 1
fi

echo "Slot jobid: ${SLOT_JOBID}"
echo
echo "=== Waiting for slot READY ==="
# Use file-existence wait (not squeue polling) — shuttle wait is the canonical path
bash "${SHUTTLE}" wait "${SLOT_JOBID}"

# 2. Get node hostname from check
NODE="$(bash "${SHUTTLE}" check "${SLOT_JOBID}" 2>&1 | awk -F= '/node=/{print $2; exit}' | tr -d ' ')"
if [[ -z "${NODE}" ]]; then
  echo "ERROR: failed to determine compute node from check output" >&2
  bash "${SHUTTLE}" check "${SLOT_JOBID}" 2>&1 | tail -5
  exit 1
fi
echo "Compute node: ${NODE}"

# 3. SSH-launch the atlas slot-runner on the compute node
REMOTE_LOG="/tmp/atlas_slot_${SLOT_JOBID}.log"
RESULTS_DIR="/rds/user/wr286/hpc-work/tidal/hpc_results/${SLOT_JOBID}"

echo
echo "=== Launching run_atlas_slot.sh ${THEORIES_STR} on ${NODE} ==="
# shellcheck disable=SC2029
ssh csd3 "ssh ${NODE} 'cd /rds/user/wr286/hpc-work/tidal && \
  export SLURM_JOB_ID=${SLOT_JOBID} && \
  export RESULTS_DIR=${RESULTS_DIR} && \
  nohup bash scripts/hpc_submit_drafts/v3_atlas/run_atlas_slot.sh ${THEORIES_STR} > ${REMOTE_LOG} 2>&1 & disown && echo launched'"

echo
echo "=== Slot ${SLOT_JOBID} running atlas for: ${THEORIES_STR} ==="
echo "Monitor with:"
echo "  ssh csd3 \"ssh ${NODE} 'tail -f ${REMOTE_LOG}'\""
echo "Wait for completion:"
echo "  until ssh csd3 \"ssh ${NODE} 'test -f ${RESULTS_DIR}/BATCH_DONE'\"; do sleep 120; done"
echo "Pull results:"
echo "  bash scripts/hpc_shuttle.sh pull ${SLOT_JOBID} --src ${RESULTS_DIR}"
echo
echo "SLOT_JOBID=${SLOT_JOBID}"
echo "NODE=${NODE}"
echo "REMOTE_LOG=${REMOTE_LOG}"
echo "RESULTS_DIR=${RESULTS_DIR}"
