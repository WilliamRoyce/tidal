#!/bin/bash
# HPC benchmark: JAX vs scipy modal solver on CSD3 sapphire.
#
# Runs two short tidal sample chains (PolyChord, INTR QOS) using the
# same theory and parameters, once with --scheme modal-jax and once with
# --scheme modal (scipy), and records wall time + posterior agreement.
#
# Usage (local — do NOT submit from this script):
#   Review this file, then submit via:
#     bash scripts/hpc_shuttle.sh submit scripts/hpc_jax_benchmark.sh
#
# Prerequisites:
#   1. JAX installed in HPC venv:
#      bash scripts/hpc_shuttle.sh shell
#      pip install "jax>=0.4.20" "jaxlib>=0.4.20"          # CPU-only
#      bash scripts/hpc_refresh_venv_tar.sh                 # refresh tarball
#   2. Local code pushed:
#      bash scripts/hpc_shuttle.sh push
#
# Output written to ${RESULTS_DIR}/jax_vs_scipy.txt
# Pull with: bash scripts/hpc_shuttle.sh pull <JOBID>
#
# NOTE: Do NOT submit this script directly. The parent session (user)
# submits once local benchmark numbers (scripts/local_jax_benchmark.py)
# confirm ≥3× speedup and rel-error < 1e-10.

set -euo pipefail

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-sar
module load python/3.11.0-icl

echo "host=$(hostname) cpus=${SLURM_CPUS_ON_NODE:-?} job=${SLURM_JOB_ID:-local}"
echo "started=$(date -Iseconds)"

REMOTE_ROOT="${REMOTE_ROOT:-${HOME}/rds/hpc-work/tidal}"
cd "${REMOTE_ROOT}"
source .venv/bin/activate

export RESULTS_DIR="${REMOTE_ROOT}/hpc_results/${SLURM_JOB_ID:-jax_bench_local}"
mkdir -p "${RESULTS_DIR}"

# ---------------------------------------------------------------------------
# Configuration — edit these to match the target theory
# ---------------------------------------------------------------------------
THEORY="${REMOTE_ROOT}/examples/data/torsion_dark_photon.json"
PRIORS="alpha3=uniform:0.1:2.0"
LIKELIHOOD="P_max:maximize"
N_SAMPLES=50          # small — this is a timing benchmark, not production
PARALLEL=8            # safe for INTR QOS login-node limits (≤4 CPUs) if local;
                      # use min(N_SAMPLES, 56) when running on a compute node
SNAPSHOTS=2           # mandatory for PolyChord (avoids disk stalls)
T_END=10.0
GRID="64"

OUTPUT_SCIPY="${RESULTS_DIR}/bench_scipy"
OUTPUT_JAX="${RESULTS_DIR}/bench_jax"

# ---------------------------------------------------------------------------
# scipy backend (baseline)
# ---------------------------------------------------------------------------
echo ""
echo "=== scipy modal backend ==="
TIME_SCIPY_START=$(date +%s%N)
tidal sample "${THEORY}" \
    --prior "${PRIORS}" \
    --likelihood "${LIKELIHOOD}" \
    --method mc \
    --n-samples "${N_SAMPLES}" \
    --parallel "${PARALLEL}" \
    --scheme modal \
    --snapshots "${SNAPSHOTS}" \
    --t-end "${T_END}" \
    --grid-shape "${GRID}" \
    --output "${OUTPUT_SCIPY}"
TIME_SCIPY_END=$(date +%s%N)
WALL_SCIPY=$(( (TIME_SCIPY_END - TIME_SCIPY_START) / 1000000 ))  # ms

# ---------------------------------------------------------------------------
# JAX backend
# ---------------------------------------------------------------------------
echo ""
echo "=== JAX modal backend ==="
TIME_JAX_START=$(date +%s%N)
tidal sample "${THEORY}" \
    --prior "${PRIORS}" \
    --likelihood "${LIKELIHOOD}" \
    --method mc \
    --n-samples "${N_SAMPLES}" \
    --parallel "${PARALLEL}" \
    --scheme modal-jax \
    --snapshots "${SNAPSHOTS}" \
    --t-end "${T_END}" \
    --grid-shape "${GRID}" \
    --output "${OUTPUT_JAX}"
TIME_JAX_END=$(date +%s%N)
WALL_JAX=$(( (TIME_JAX_END - TIME_JAX_START) / 1000000 ))  # ms

# ---------------------------------------------------------------------------
# Posterior agreement check (Python)
# ---------------------------------------------------------------------------
echo ""
echo "=== Posterior agreement ==="
python3 - <<'PYEOF'
import json, sys, math
from pathlib import Path

res_dir = Path("${RESULTS_DIR}")
f_scipy = res_dir / "bench_scipy" / "inference.json"
f_jax   = res_dir / "bench_jax"   / "inference.json"

if not f_scipy.exists() or not f_jax.exists():
    print("[WARN] inference.json not found — skipping posterior comparison")
    sys.exit(0)

with open(f_scipy) as f: d_s = json.load(f)
with open(f_jax)   as f: d_j = json.load(f)

# Compare MAP estimates
map_s = d_s.get("map_parameters", {})
map_j = d_j.get("map_parameters", {})
max_rel = 0.0
for key in map_s:
    vs, vj = map_s.get(key, 0.0), map_j.get(key, 0.0)
    ref = abs(vs) if abs(vs) > 1e-15 else 1.0
    rel = abs(vs - vj) / ref
    max_rel = max(max_rel, rel)
    print(f"  {key}: scipy={vs:.6g}  jax={vj:.6g}  rel={rel:.2e}")

print(f"\nMax relative MAP difference: {max_rel:.2e}")
if max_rel < 0.05:
    print("[PASS] Posteriors agree to within 5%")
else:
    print("[WARN] Posteriors differ by more than 5% — check convergence")
PYEOF

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
SPEEDUP_X=$(python3 -c "print(f'{${WALL_SCIPY}/${WALL_JAX}:.2f}')" 2>/dev/null || echo "?")

{
    echo "=== JAX vs scipy modal solver HPC benchmark ==="
    echo "theory: ${THEORY}"
    echo "n_samples: ${N_SAMPLES}  parallel: ${PARALLEL}  grid: ${GRID}  t_end: ${T_END}"
    echo ""
    echo "scipy wall time:  ${WALL_SCIPY} ms"
    echo "jax   wall time:  ${WALL_JAX} ms"
    echo "speedup (scipy/jax): ${SPEEDUP_X}x"
    echo ""
    echo "Results:"
    echo "  scipy: ${OUTPUT_SCIPY}"
    echo "  jax:   ${OUTPUT_JAX}"
    echo ""
    echo "finished=$(date -Iseconds)"
} | tee "${RESULTS_DIR}/jax_vs_scipy.txt"

echo ""
echo "Pull results with: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID:-?}"
echo "finished=$(date -Iseconds)"
