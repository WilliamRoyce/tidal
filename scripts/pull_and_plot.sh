#!/usr/bin/env bash
# pull_and_plot.sh — atomically pull an HPC nested-sampling chain AND generate
# its corner plot.  Mandatory post-job step per docs/V3_PHASE_D_DESIGN.md
# §"Chain-pull workflow" — a chain without a corner plot is half-pulled.
#
# Usage:
#   bash scripts/pull_and_plot.sh <jobid> [--src REMOTE_PATH] [--label LABEL]
#                                  [--v2 V2_JOBID]
#
# What it does:
#   1. hpc_shuttle.sh pull <jobid> [--src ...]
#   2. uv run tidal plot hpc_results/<jobid> --type corner --show-rejected-prior
#      --output hpc_results/<jobid>/corner_v3.png
#   3. Print a one-line summary (log Z, ESS, MAP) from inference.json
#   4. Print run_status breakdown from results.csv
#   5. If --v2 given: generate docs/comparison/<label>_v2_v3.md
#
# Exit non-zero if any step fails — the chain is not "delivered" until all
# three deliverables exist locally.
#
# Examples:
#   bash scripts/pull_and_plot.sh 29149987
#   bash scripts/pull_and_plot.sh 29189966 --label "stage_a_amp" --v2 28474676
#   bash scripts/pull_and_plot.sh 29189748 --src /rds/user/wr286/hpc-work/tidal/hpc_results/29189748/d1_amp_v3

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
    sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
    exit "${1:-0}"
}

jobid=""
src=""
label=""
v2_jobid=""
no_pull=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage ;;
        --src)     src="$2"; shift 2 ;;
        --label)   label="$2"; shift 2 ;;
        --v2)      v2_jobid="$2"; shift 2 ;;
        --no-pull) no_pull="1"; shift ;;
        *)         jobid="$1"; shift ;;
    esac
done
[[ -n "$jobid" ]] || { echo "error: jobid required" >&2; usage 1; }

cd "$REPO_ROOT"

# Step 1: pull (skip if --no-pull and chain already local)
if [[ -n "$no_pull" ]]; then
    echo "==> [1/4] --no-pull: assuming hpc_results/${jobid}/ already exists"
else
    echo "==> [1/4] Pulling chain artifacts for $jobid"
    if [[ -n "$src" ]]; then
        bash scripts/hpc_shuttle.sh pull "$jobid" --src "$src"
    else
        bash scripts/hpc_shuttle.sh pull "$jobid"
    fi
fi

dst="hpc_results/${jobid}"
[[ -f "${dst}/inference.json" ]] || {
    echo "error: ${dst}/inference.json missing after pull; the chain didn't complete or --src is wrong" >&2
    exit 2
}

# Step 2: corner plot
echo "==> [2/4] Generating corner plot"
corner_png="${dst}/corner_v3.png"
uv run tidal plot "$dst" --type corner --show-rejected-prior --output "$corner_png"
[[ -f "$corner_png" ]] || {
    echo "error: corner plot was not written to $corner_png" >&2
    exit 3
}

# Step 3: chain-health summary
echo "==> [3/4] Chain summary"
uv run python -c "
import json
import csv
from collections import Counter
from pathlib import Path

d = json.load(open('${dst}/inference.json'))
print(f'  log Z       : {d[\"log_evidence\"]:+.3f} ± {d[\"log_evidence_err\"]:.3f}')
print(f'  ESS         : {d[\"effective_sample_size\"]:.0f} / {d[\"n_samples\"]} samples')
print(f'  Joint D_KL  : {d[\"parameter_importance\"][\"d_kl\"]:.2f} nats')
print(f'  MAP         : {d[\"map_estimate\"]}')
print(f'  Per-coupling D_KL: {d[\"parameter_importance\"][\"marginal_d_kl\"]}')
csv_path = Path('${dst}/results.csv')
if csv_path.exists():
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    if rows and 'run_status' in rows[0]:
        counts = Counter(r['run_status'] for r in rows)
        print(f'  run_status breakdown:')
        for k, v in counts.most_common():
            pct = 100*v/len(rows)
            flag = '  ⚠' if k != 'success' and k != 'below_noise_floor' and pct > 5 else ''
            print(f'    {k:25s}: {v:5d} ({pct:.1f}%){flag}')
"

# Step 4 (optional): v2 vs v3 comparison
if [[ -n "$v2_jobid" ]]; then
    echo "==> [4/4] Generating v2 vs v3 comparison vs job $v2_jobid"
    [[ -f "hpc_results/${v2_jobid}/inference.json" ]] || {
        echo "error: hpc_results/${v2_jobid}/inference.json missing; pull the v2 reference first" >&2
        exit 4
    }
    comp_label="${label:-${jobid}_vs_${v2_jobid}}"
    comp_out="docs/comparison/${comp_label}_v2_v3.md"
    uv run python scripts/v3_v2_comparison.py \
        --v2 "hpc_results/${v2_jobid}" \
        --v3 "$dst" \
        --label "$comp_label" \
        --output "$comp_out"
else
    echo "==> [4/4] (no --v2 given; skipping comparison)"
fi

echo
echo "==> Done: $jobid"
echo "  Chain dir   : $dst/"
echo "  Corner plot : $corner_png"
[[ -n "$v2_jobid" ]] && echo "  Comparison  : $comp_out"
