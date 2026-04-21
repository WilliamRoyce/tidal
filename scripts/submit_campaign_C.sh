#!/bin/bash
# Submit Stage C: R²-PGT perturbative b5 (T3), both perturbative orders.
#
# Gate: Stage C gate must be cleared before running (see CAMPAIGN.md).
# Note: b5 expected null in h_5→a_1 channel (see issue #299, #300).
# These runs confirm D_KL(b5,ord=1) vs D_KL(b5,ord=0) — structural comparison.
#
# Usage:
#   bash scripts/submit_campaign_C.sh

set -euo pipefail
TEMPLATE="scripts/hpc_templates/polychord_intr.sbatch"

echo "=== Stage C: R²-PGT b5 perturbative (T3) — both orders ==="

for ORDER in 0 1; do
  bash scripts/hpc_shuttle.sh submit \
    --template "$TEMPLATE" \
    --nodes 1 --ntasks 16 --time 01:00:00 --name "b5_ord${ORDER}" \
    --cmd "OUTPUT_DIR=\${REMOTE_ROOT}/results/b5_ord${ORDER}_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein.json \
  --param kappa=1.0 --param B0=0.01 \
  --perturbative-order ${ORDER} \
  --prior 'b5=log_uniform:1e-4:1e-2' \
  --prior 'alpha1=uniform:-1:1' \
  --prior 'alpha2=uniform:-1:1' \
  --prior 'alpha3=log_uniform:0.05:2' \
  --likelihood 'P_max:maximize' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 400 \
  --num-repeats 5 --precision-criterion 0.01 --no-clustering \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_ord${ORDER}.png"
done

echo ""
echo "=== Stage C: 2 jobs submitted (order=0 and order=1). ==="
echo "Success criterion: D_KL(b5,ord=1) − D_KL(b5,ord=0) > 0.01 nats"
echo "Expected: both ~ 0 (structural null — see issue #299)"
