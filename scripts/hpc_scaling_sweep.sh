#!/usr/bin/env bash
# HPC scaling test: run the plasma-Gertsenshtein C₀ sweep at multiple
# --parallel values within a single sbatch allocation, printing wall
# time for each. Produces a scaling curve for HPC parallelism diagnosis.
#
# Designed to be invoked from the sbatch --cmd substitution, so it
# takes no arguments and uses repo-relative paths.
set -euo pipefail

OUT_ROOT="/rds/user/wr286/hpc-work/tidal/hpc_smoke"

for P in 1 8 16 32 64 90; do
    echo ""
    echo "=== SCALING parallel=${P} ==="
    START=$(date +%s.%N)

    tidal sweep examples/data/gertsenshtein_proca.json \
        --sweep "B0=1e-4,3e-4,1e-3" \
        --sweep "mA2=0.0:0.1:30" \
        --param kappa=1.0 \
        --measure conversion,conservation \
        --source h_7 --target a_2 \
        --grid-shape 256 --bounds 0:100 --periodic \
        --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 \
        --ic-amplitude 0.1 \
        --t-end 30 \
        --parallel "${P}" --force \
        --output "${OUT_ROOT}/scaling_p${P}" 2>&1 | tail -5

    END=$(date +%s.%N)
    ELAPSED=$(echo "${END} - ${START}" | bc)
    printf "\nSCALING_RESULT parallel=%s wall=%.2f s\n" "${P}" "${ELAPSED}"
done

echo ""
echo "=== SCALING TEST COMPLETE ==="
