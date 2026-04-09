#!/usr/bin/env bash
# Phase 1: 1D Parameter Slices for General Nonminimal PGT Model
#
# Sweeps each of the 9 torsion/coupling parameters independently
# while holding all others at the baseline.
#
# Baseline: beta1=0, beta2=-0.6, beta3=0.5, xi=1, delta1=0, chi=0, zeta1=0, zeta2=0, zeta3=0
# Fixed: kappa=1, B0=0.01, amplitude=0.001 (amp/B0 = 0.1)
#
# Each sweep: 50 arctan-mapped LHS samples spanning (-infty, +infty)
# Measurement: peak_conversion (h_5 -> a_1)
# Tachyonic pre-check: skips unstable points automatically
#
# Usage: bash sweep_phase1_1d.sh [OUTPUT_DIR]

set -euo pipefail

SPEC="examples/data/torsion_gertsenshtein_general_nonminimal.json"
OUT="${1:-/tmp/tidal_sweeps/phase1}"

COMMON_ARGS=(
    --measure peak_conversion
    --source h_5 --target a_1
    --grid-shape 256 --bounds 0:100 --periodic
    --t-end 50 --snapshots 50
    --ic plane-wave --ic-amplitude 1e-3 --ic-component h_5
    --param kappa=1 --param B0=0.01
    --parallel 4
    --sweep-strategy latin_hypercube --n-samples 50
)

BASELINE_PARAMS=(
    --param beta1=0 --param beta2=-0.6 --param beta3=0.5 --param xi=1.0
    --param delta1=0 --param chi=0 --param zeta1=0 --param zeta2=0 --param zeta3=0
)

PARAMS=(beta1 beta2 beta3 xi delta1 chi zeta1 zeta2 zeta3)

for P in "${PARAMS[@]}"; do
    echo "=== Sweeping $P ==="
    uv run tidal sweep "$SPEC" \
        --sweep "${P}=50:arctan" \
        "${COMMON_ARGS[@]}" \
        "${BASELINE_PARAMS[@]}" \
        --output "${OUT}/${P}" \
        --resume
    echo ""
done

echo "Phase 1 complete. Results in: $OUT"
