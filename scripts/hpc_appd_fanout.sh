#!/usr/bin/env bash
# Fan-out the six App D calibration benchmarks in parallel on a single
# HPC compute node, writing canonical JSON to
# benchmark_results/canonical/. Invoked from an interactive INTR slot
# attached via scripts/hpc_shuttle.sh attach (or directly via ssh from
# the devcontainer).
#
# Each benchmark already manages its own internal parallelism for tidal
# sweep where applicable. Total wall is bounded by the slowest single
# benchmark (likely boccaletti_calibration with its 40-point B0 sweep
# plus the multi-N convergence). All six run concurrently with `&` and
# the script blocks on `wait`.

set -uo pipefail

REPO_ROOT="${REPO_ROOT:-/rds/user/wr286/hpc-work/tidal}"
cd "$REPO_ROOT"

# Ensure venv is active (idempotent if already activated by attach).
if [[ -z "${VIRTUAL_ENV:-}" ]]; then
  source .venv/bin/activate
fi

LOG_DIR="${LOG_DIR:-${REPO_ROOT}/benchmark_results/canonical/_logs}"
mkdir -p "$LOG_DIR" "${REPO_ROOT}/benchmark_results/canonical"

echo "=== App D calibration fan-out ==="
echo "started: $(date -Is)"
echo "node:    $(hostname)"
echo "logs:    $LOG_DIR"

# Six benchmarks in parallel. Per-script parallelism tuned to the
# internal work each does. Sapphire has 112 cores.
python scripts/benchmarks/boccaletti_calibration.py --parallel 16 > "$LOG_DIR/boccaletti.log" 2>&1 &
PID_BOCC=$!
python scripts/benchmarks/torsion_limit.py --parallel 8 > "$LOG_DIR/torsion.log" 2>&1 &
PID_TORS=$!
python scripts/benchmarks/cross_backend.py > "$LOG_DIR/cross_backend.log" 2>&1 &
PID_CROS=$!
python scripts/benchmarks/proca_dispersion.py > "$LOG_DIR/proca.log" 2>&1 &
PID_PROC=$!
python scripts/benchmarks/time_integration_order.py > "$LOG_DIR/tio.log" 2>&1 &
PID_TIO=$!
python scripts/benchmarks/conservation_audit.py > "$LOG_DIR/conservation.log" 2>&1 &
PID_CONS=$!

echo "PIDs: bocc=$PID_BOCC torsion=$PID_TORS cross_backend=$PID_CROS proca=$PID_PROC tio=$PID_TIO conservation=$PID_CONS"

# Wait for each, collecting exit codes.
declare -A FAIL=()
for name in BOCC TORS CROS PROC TIO CONS; do
  pid_var="PID_$name"
  pid="${!pid_var}"
  wait "$pid"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    FAIL[$name]=$rc
  fi
  echo "  [$name] pid=$pid rc=$rc"
done

echo
echo "=== completed ==="
echo "ended: $(date -Is)"

# Summary of generated canonical JSON.
ls -la "${REPO_ROOT}/benchmark_results/canonical/"*.json 2>/dev/null || echo "(no canonical JSON found)"

# Exit non-zero if any benchmark failed.
if [[ ${#FAIL[@]} -gt 0 ]]; then
  echo "FAILURES: ${!FAIL[*]}"
  exit 1
fi

echo "=== ALL OK ==="
