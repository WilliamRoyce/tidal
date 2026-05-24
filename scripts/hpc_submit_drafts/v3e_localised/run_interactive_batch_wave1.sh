#!/usr/bin/env bash
# Phase E Wave 1+2 batch runner — chains multiple theories' amp+sup in one INTR slot.
#
# Per-theory: amp+sup parallel (32 ranks each = 64 cores), ~5 min based on
# per-eval benchmark at the new compact geometry (L=100, grid=128, t_end=40).
# Five theories sequential = ~25 min, well within 60-min INTR budget.
#
# Run INSIDE an interactive sapphire allocation:
#   bash scripts/hpc_shuttle.sh interactive
#   bash scripts/hpc_shuttle.sh attach <jobid>   # interactive
# OR launch via SSH:
#   ssh csd3 "ssh <node> 'cd /rds/.../tidal && export SLURM_JOB_ID=<jobid> && \
#       export RESULTS_DIR=/rds/.../hpc_results/<jobid> && \
#       bash scripts/hpc_submit_drafts/v3e_localised/run_interactive_batch_wave1.sh \
#       > /tmp/wave1.log 2>&1 &'"

set -euo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"
SCRIPT_DIR="${TIDAL_ROOT}/scripts/hpc_submit_drafts/v3e_localised"

# Theories to run in order (cheapest physics first; T7s last since it's 18D).
THEORIES=(t4 t1 t5 np t7s)

echo "=== [$(date +%H:%M:%S)] Wave-1 batch: ${THEORIES[*]} ==="
for t in "${THEORIES[@]}"; do
  echo
  echo "=== [$(date +%H:%M:%S)] === Starting theory: ${t} ==="
  bash "${SCRIPT_DIR}/run_interactive_amp_sup.sh" "${t}" || {
    echo "WARN: theory ${t} runner failed; continuing with next theory"
  }
  echo "=== [$(date +%H:%M:%S)] === Finished theory: ${t} ==="
done
echo
echo "=== [$(date +%H:%M:%S)] BATCH DONE ==="
date
