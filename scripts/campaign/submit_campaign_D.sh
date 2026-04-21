#!/bin/bash
# Submit Stage D: nonminimal survey — D1 (T4 Ricci-EM), D2 (T5 YM-PGT + sub-theories),
# D3 (T6 YM-PGT-CP). Run after Stages A-C have been pulled and reviewed.
#
# Usage:
#   bash scripts/campaign/submit_campaign_D.sh [d1|d2|d3|all]
#   Default: all

set -euo pipefail
TEMPLATE="scripts/hpc_templates/polychord_intr.sbatch"
STAGE="${1:-all}"

# ---- Stage D1: Ricci-EM / T4 (4 params) ----
submit_d1() {
  echo "=== Stage D1: Ricci-EM (T4) — paired amplify + suppress ==="
  for MODE in amplify suppress; do
    LHOOD=$( [[ $MODE == amplify ]] && echo "P_max:maximize" || echo "P_max:minimize" )
    bash scripts/hpc_shuttle.sh submit \
      --template "$TEMPLATE" \
      --nodes 1 --ntasks 16 --time 01:00:00 --name "ricciem_${MODE}" \
      --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/ricciem_${MODE}_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
  --param kappa=1.0 --param B0=0.01 \
  --prior 'alpha1=uniform:-2:2' \
  --prior 'alpha2=uniform:-2:2' \
  --prior 'alpha3=log_uniform:0.01:4' \
  --prior 'delta1=uniform:-0.5:0.5' \
  --likelihood '${LHOOD}' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 400 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_${MODE}.png"
  done
}

# ---- Stage D2: YM-PGT / T5 (11 params) + paper sub-theories ----
submit_d2() {
  echo "=== Stage D2: YM-PGT (T5) — full + sub-theories, paired ==="
  # Full YM-PGT: all 11 params free
  for MODE in amplify suppress; do
    LHOOD=$( [[ $MODE == amplify ]] && echo "P_max:maximize" || echo "P_max:minimize" )
    bash scripts/hpc_shuttle.sh submit \
      --template "$TEMPLATE" \
      --nodes 1 --ntasks 16 --time 01:00:00 --name "ympgt_${MODE}" \
      --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/ympgt_${MODE}_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
  --param kappa=1.0 --param B0=0.01 \
  --prior 'beta1=uniform:-2:2' \
  --prior 'beta2=uniform:-2:2' \
  --prior 'beta3=uniform:-2:2' \
  --prior 'xi=log_uniform:0.001:10.0' \
  --prior 'delta1=uniform:-0.5:0.5' \
  --prior 'chi=uniform:-0.5:0.5' \
  --prior 'zeta1=uniform:-0.5:0.5' \
  --prior 'zeta2=uniform:-0.5:0.5' \
  --prior 'zeta3=uniform:-0.5:0.5' \
  --likelihood '${LHOOD}' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 500 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_${MODE}.png"
  done

  # Barker-PGT sub-theory: beta1-3, xi, chi free; delta1=zeta1=zeta2=zeta3=0
  bash scripts/hpc_shuttle.sh submit \
    --template "$TEMPLATE" \
    --nodes 1 --ntasks 16 --time 01:00:00 --name "barker_amplify" \
    --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/barker_amplify_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
  --param kappa=1.0 --param B0=0.01 \
  --param delta1=0 --param zeta1=0 --param zeta2=0 --param zeta3=0 \
  --prior 'beta1=uniform:-2:2' \
  --prior 'beta2=uniform:-2:2' \
  --prior 'beta3=uniform:-2:2' \
  --prior 'xi=log_uniform:0.001:10.0' \
  --prior 'chi=uniform:-0.5:0.5' \
  --likelihood 'P_max:maximize' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 400 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_amplify.png"

  # Shapiro-PGT sub-theory: beta1-3, zeta1-3 free; xi=delta1=chi=0
  bash scripts/hpc_shuttle.sh submit \
    --template "$TEMPLATE" \
    --nodes 1 --ntasks 16 --time 01:00:00 --name "shapiro_amplify" \
    --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/shapiro_amplify_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
  --param kappa=1.0 --param B0=0.01 \
  --param xi=0 --param delta1=0 --param chi=0 \
  --prior 'beta1=uniform:-2:2' \
  --prior 'beta2=uniform:-2:2' \
  --prior 'beta3=uniform:-2:2' \
  --prior 'zeta1=uniform:-0.5:0.5' \
  --prior 'zeta2=uniform:-0.5:0.5' \
  --prior 'zeta3=uniform:-0.5:0.5' \
  --likelihood 'P_max:maximize' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 400 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_amplify.png"
}

# ---- Stage D3: YM-PGT-CP / T6 (24 params) ----
submit_d3() {
  echo "=== Stage D3: YM-PGT-CP (T6) — paired amplify + suppress ==="
  for MODE in amplify suppress; do
    LHOOD=$( [[ $MODE == amplify ]] && echo "P_max:maximize" || echo "P_max:minimize" )
    bash scripts/hpc_shuttle.sh submit \
      --template "$TEMPLATE" \
      --nodes 1 --ntasks 16 --time 03:00:00 --name "ympgtcp_${MODE}" \
      --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/ympgtcp_${MODE}_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein_parity_odd.json \
  --param kappa=1.0 --param B0=0.01 \
  --prior 'beta1=uniform:-2:2' \
  --prior 'beta2=uniform:-2:2' \
  --prior 'beta3=uniform:-2:2' \
  --prior 'xi=log_uniform:0.001:10.0' \
  --prior 'delta1=uniform:-0.5:0.5' \
  --prior 'chi=uniform:-0.5:0.5' \
  --prior 'zeta1=uniform:-0.5:0.5' \
  --prior 'zeta2=uniform:-0.5:0.5' \
  --prior 'zeta3=uniform:-0.5:0.5' \
  --prior 'd14=uniform:-0.5:0.5' \
  --prior 'd15=uniform:-0.5:0.5' \
  --prior 'd16=uniform:-0.5:0.5' \
  --prior 'd17=uniform:-0.5:0.5' \
  --prior 'd19=uniform:-0.5:0.5' \
  --prior 'd20=uniform:-0.5:0.5' \
  --prior 'd21=uniform:-0.5:0.5' \
  --prior 'zt1=uniform:-0.5:0.5' \
  --prior 'zt2=uniform:-0.5:0.5' \
  --prior 'zt3=uniform:-0.5:0.5' \
  --prior 'zt4=uniform:-0.5:0.5' \
  --prior 'zt5=uniform:-0.5:0.5' \
  --prior 'zt6=uniform:-0.5:0.5' \
  --likelihood '${LHOOD}' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 500 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_${MODE}.png"
  done
}

case "$STAGE" in
  d1)  submit_d1 ;;
  d2)  submit_d2 ;;
  d3)  submit_d3 ;;
  all) submit_d1; submit_d2; submit_d3 ;;
  *)   echo "Usage: $0 [d1|d2|d3|all]"; exit 1 ;;
esac

echo ""
echo "=== Stage D jobs submitted. Record job IDs in CAMPAIGN.md. ==="
