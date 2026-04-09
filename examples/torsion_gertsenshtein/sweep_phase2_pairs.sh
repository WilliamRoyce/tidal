#!/usr/bin/env bash
# Phase 2: Pairwise 2D Parameter Interactions
#
# 10 priority pairs of parameters swept simultaneously.
# 200 arctan-mapped LHS samples per pair.
#
# Pairs: coupling x coupling, coupling x mass, coupling x kinetics
#
# Usage: bash sweep_phase2_pairs.sh [OUTPUT_DIR]

set -euo pipefail

SPEC="examples/data/torsion_gertsenshtein_general_nonminimal.json"
OUT="${1:-/tmp/tidal_sweeps/phase2}"

COMMON_ARGS=(
    --measure peak_conversion
    --source h_5 --target a_1
    --grid-shape 256 --bounds 0:100 --periodic
    --t-end 50 --snapshots 50
    --ic plane-wave --ic-amplitude 1e-3 --ic-component h_5
    --param kappa=1 --param B0=0.01
    --parallel 4
    --sweep-strategy latin_hypercube --n-samples 200
)

BASELINE_PARAMS=(
    --param beta1=0 --param beta2=-0.6 --param beta3=0.5 --param xi=1.0
    --param delta1=0 --param chi=0 --param zeta1=0 --param zeta2=0 --param zeta3=0
)

# 10 priority pairs
PAIRS=(
    "delta1 zeta1"
    "delta1 zeta2"
    "delta1 zeta3"
    "zeta1 zeta2"
    "zeta1 zeta3"
    "zeta2 zeta3"
    "delta1 beta2"
    "zeta1 beta2"
    "delta1 xi"
    "zeta1 xi"
)

for PAIR in "${PAIRS[@]}"; do
    P1=$(echo "$PAIR" | cut -d' ' -f1)
    P2=$(echo "$PAIR" | cut -d' ' -f2)
    echo "=== Sweeping $P1 x $P2 ==="
    uv run tidal sweep "$SPEC" \
        --sweep "${P1}=20:arctan" --sweep "${P2}=20:arctan" \
        "${COMMON_ARGS[@]}" \
        "${BASELINE_PARAMS[@]}" \
        --output "${OUT}/${P1}_${P2}" \
        --resume
    echo ""
done

echo "Phase 2 complete. Results in: $OUT"
