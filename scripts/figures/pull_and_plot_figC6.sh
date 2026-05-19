#!/usr/bin/env bash
# Pull the (possibly-in-progress) jac_speedup.json from CSD3 and regenerate
# the figC6 PDF.  Safe to run while the HPC benchmark job is still running —
# the script writes per-config checkpoints to the same path on every config
# completion, so each pull captures a strictly-larger snapshot.
#
# Usage:  bash scripts/figures/pull_and_plot_figC6.sh
#         make figC6-pull
#
# Outputs:
#   benchmark_results/canonical/jac_speedup.json   (snapshot, overwritten)
#   manuscript/figures/figC6_jac_speedup.pdf       (regenerated)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

REMOTE_PATH="csd3:/rds/user/wr286/hpc-work/tidal/benchmark_results/canonical/jac_speedup.json"
LOCAL_PATH="benchmark_results/canonical/jac_speedup.json"
SUMMARY="/tmp/figC6_summary.txt"

echo "==> pulling $REMOTE_PATH"
scp -q "$REMOTE_PATH" "$LOCAL_PATH"

echo "==> snapshot summary:"
uv run python - <<'PY' | tee "$SUMMARY"
import json, math
from pathlib import Path

p = Path("benchmark_results/canonical/jac_speedup.json")
d = json.loads(p.read_text())
rows = d.get("results", [])
print(f"  {len(rows)} configs in snapshot")
by_theory: dict[str, int] = {}
for r in rows:
    by_theory[r["theory"]] = by_theory.get(r["theory"], 0) + 1
for t, n in sorted(by_theory.items()):
    print(f"    {t:<22} {n:>3} configs")

# Speedup summary in the dense tier
print("  speedup_auto_vs_fd by tier:")
for tier in ("dense", "sparse", "gmres"):
    vals = [
        r["speedup_auto_vs_fd"] for r in rows
        if r.get("auto_tier") == tier
        and isinstance(r.get("speedup_auto_vs_fd"), (int, float))
        and not math.isnan(r["speedup_auto_vs_fd"])
    ]
    if vals:
        print(f"    {tier:<6}: n={len(vals)}  min={min(vals):.2f}x  max={max(vals):.2f}x  median={sorted(vals)[len(vals)//2]:.2f}x")
PY

echo "==> regenerating manuscript/figures/figC6_jac_speedup.pdf"
uv run python scripts/figures/figC6_jac_speedup.py

echo "==> done"
