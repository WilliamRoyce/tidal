#!/usr/bin/env bash
# Phase 3: Full 9D Latin Hypercube Survey
#
# All 9 torsion/coupling parameters swept simultaneously.
# 1000 arctan-mapped LHS samples in 9D space.
#
# Usage: bash sweep_phase3_9d.sh [OUTPUT_DIR]

set -euo pipefail

SPEC="examples/data/torsion_gertsenshtein_general_nonminimal.json"
OUT="${1:-/tmp/tidal_sweeps/phase3}"

uv run tidal sweep "$SPEC" \
    --sweep "beta1=20:arctan" --sweep "beta2=20:arctan" --sweep "beta3=20:arctan" \
    --sweep "xi=20:arctan" \
    --sweep "delta1=20:arctan" --sweep "chi=20:arctan" \
    --sweep "zeta1=20:arctan" --sweep "zeta2=20:arctan" --sweep "zeta3=20:arctan" \
    --measure peak_conversion \
    --source h_5 --target a_1 \
    --grid-shape 256 --bounds 0:100 --periodic \
    --t-end 50 --snapshots 50 \
    --ic plane-wave --ic-amplitude 1e-3 --ic-component h_5 \
    --param kappa=1 --param B0=0.01 \
    --sweep-strategy latin_hypercube --n-samples 1000 \
    --parallel 4 \
    --output "${OUT}/lhs1000" \
    --resume

echo "Phase 3 complete. Results in: $OUT"
