#!/usr/bin/env bash
# ==============================================================================
# run_lane.sh -- run one R-C Wolfram step on the single-kernel lane
# ==============================================================================
# Usage: run_lane.sh <step> [key=value ...]
# Steps (script, timeout): probe_load (15m) · a2 (30m) · a1 (45m) · b (45m) ·
#   c (30m) · d (90m) · dl (30m) · e (45m) · e2 (15m) · t1 (20m) · f2 (20m) ·
#   f3 (15m) · f5 (15m) · f6 (15m) · f7 (45m) · f8 (20m)
# Discipline: refuses if a kernel is live (the lane hook cannot see a bare .wls);
# QT_QPA_PLATFORM=offscreen; runs from a throwaway cwd; hard timeout; the
# transcript is scrubbed of the home directory and repo path at print time and
# saved under third_party/perturbations_runs/<step>/<utc>/; success is read from
# the RC_<STEP>_DONE sentinel, never from the exit status (#561).
# ==============================================================================
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="$REPO_ROOT/scripts/research/perturbations"
STEP="${1:-}"; shift || true
declare -A SCRIPT=( [probe_load]=probe_load.wls [a2]=repro_a2_tensor_eom.wls [a1]=repro_a1_tensor_action.wls
  [b]=repro_b_mb_scalars.wls [c]=probe_c_signature.wls [d]=probe_d_torsion.wls [dl]=probe_d_limits.wls [e]=probe_e_xmag.wls [e2]=probe_e2_xmag_induced.wls
  [t1]=tier1_xmag.wls [f2]=f2_hazards.wls [f3]=f3_changecurvature.wls [f5]=f5_sign_audit.wls [f6]=f6_psalter_signs.wls
  [f7]=f7_epsilon_and_import_map.wls [f8]=f8_contractmetric.wls )
declare -A TMO=( [probe_load]=15m [a2]=30m [a1]=45m [b]=45m [c]=30m [d]=90m [dl]=30m [e]=45m [e2]=15m [t1]=20m [f2]=20m [f3]=15m [f5]=15m [f6]=15m [f7]=45m [f8]=20m )
[[ -n "$STEP" && -n "${SCRIPT[$STEP]:-}" ]] || { echo "usage: $0 <${!SCRIPT[*]}> [key=value ...]"; exit 2; }
if pgrep -x wolframscript >/dev/null || pgrep -x WolframKernel >/dev/null; then
    echo "RC_LANE_BUSY=a Wolfram kernel is live; refusing (single license)"; exit 3
fi
UTC="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${RUN_DIR:-$REPO_ROOT/third_party/perturbations_runs/$STEP/$UTC}"
mkdir -p "$RUN_DIR"
scrub() { sed -e "s|$HOME|~|g" -e "s|$REPO_ROOT|<repo>|g"; }
WORK="$(mktemp -d)"
START=$(date +%s)
( cd "$WORK" && QT_QPA_PLATFORM=offscreen timeout --signal=INT --kill-after=60 "${TMO[$STEP]}" \
    wolframscript -file "$HERE/wolfram/${SCRIPT[$STEP]}" "outdir=$RUN_DIR" "$@" ) 2>&1 | scrub | tee "$RUN_DIR/transcript.txt"
STATUS=${PIPESTATUS[0]}
END=$(date +%s)
rm -rf "$WORK"
echo "RC_${STEP^^}_WALL_S=$((END-START))" | tee -a "$RUN_DIR/transcript.txt"
echo "RC_${STEP^^}_EXIT_STATUS=$STATUS (informational; the sentinel is the verdict)" | tee -a "$RUN_DIR/transcript.txt"
if grep -q "^RC_${STEP^^}_DONE" "$RUN_DIR/transcript.txt"; then
    echo "RC_${STEP^^}_SENTINEL=present"; echo "RC_RUN_DIR=$(echo "$RUN_DIR" | scrub)"
else
    echo "RC_${STEP^^}_SENTINEL=MISSING (step did not complete)"; echo "RC_RUN_DIR=$(echo "$RUN_DIR" | scrub)"; exit 1
fi
