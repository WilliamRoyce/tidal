#!/usr/bin/env bash
# Fan-out the eight App-D calibration benchmarks in parallel on a single
# HPC compute node, writing canonical JSON to
# benchmark_results/canonical/. Invoked from an interactive INTR slot
# attached via scripts/hpc_shuttle.sh attach (or directly via ssh from
# the devcontainer).
#
# Each benchmark already manages its own internal parallelism for tidal
# sweep where applicable. Total wall is bounded by the slowest single
# benchmark. All eight run concurrently with `&` and the script blocks
# on `wait`.

set -uo pipefail

REPO_ROOT="${REPO_ROOT:-/rds/user/wr286/hpc-work/tidal}"
cd "$REPO_ROOT"

if [[ -z "${VIRTUAL_ENV:-}" ]]; then
  source .venv/bin/activate
fi

LOG_DIR="${LOG_DIR:-${REPO_ROOT}/benchmark_results/canonical/_logs}"
mkdir -p "$LOG_DIR" "${REPO_ROOT}/benchmark_results/canonical"

echo "=== App-D calibration fan-out (v2) ==="
echo "started: $(date -Is)"
echo "node:    $(hostname)"
echo "logs:    $LOG_DIR"

# Eight benchmarks in parallel. Sapphire has 112 cores.
python scripts/benchmarks/boccaletti_uniform_rich.py --parallel 24 > "$LOG_DIR/boccaletti_uniform.log" 2>&1 &
PID_BUR=$!
python scripts/benchmarks/boccaletti_localised.py > "$LOG_DIR/boccaletti_localised.log" 2>&1 &
PID_BLOC=$!
python scripts/benchmarks/forced_oscillator.py > "$LOG_DIR/forced_oscillator.log" 2>&1 &
PID_FO=$!
python scripts/benchmarks/parker_simon_flrw.py > "$LOG_DIR/parker_simon.log" 2>&1 &
PID_PS=$!
python scripts/benchmarks/fv_cdt_equivalence.py > "$LOG_DIR/fv_cdt.log" 2>&1 &
PID_FVCDT=$!
python scripts/benchmarks/multi_method.py > "$LOG_DIR/multi_method.log" 2>&1 &
PID_MM=$!
python scripts/benchmarks/convergence_rich.py > "$LOG_DIR/convergence_rich.log" 2>&1 &
PID_CR=$!
python scripts/benchmarks/conservation_audit_full.py > "$LOG_DIR/conservation.log" 2>&1 &
PID_CONS=$!

echo "PIDs:"
echo "  boccaletti_uniform=$PID_BUR"
echo "  boccaletti_localised=$PID_BLOC"
echo "  forced_oscillator=$PID_FO"
echo "  parker_simon=$PID_PS"
echo "  fv_cdt=$PID_FVCDT"
echo "  multi_method=$PID_MM"
echo "  convergence_rich=$PID_CR"
echo "  conservation=$PID_CONS"

declare -A FAIL=()
for name in BUR BLOC FO PS FVCDT MM CR CONS; do
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
ls -la "${REPO_ROOT}/benchmark_results/canonical/"*.json 2>/dev/null || echo "(no canonical JSON found)"

if [[ ${#FAIL[@]} -gt 0 ]]; then
  echo "FAILURES: ${!FAIL[*]}"
  exit 1
fi

echo "=== ALL OK ==="
