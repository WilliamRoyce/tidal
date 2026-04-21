#!/bin/bash
# Submit Stage A (T1 Dark-Photon-Plasma, paired amplify+suppress)
# and Stage B (T2 Einstein-Cartan null) simultaneously.
#
# Prerequisites (run once if not already done):
#   bash scripts/hpc_shuttle.sh push
#
# Usage:
#   bash scripts/campaign/submit_campaign_AB.sh
#
# After completion, record the job IDs printed in CAMPAIGN.md HPC Job Log.

set -euo pipefail
# INTR QOS has MaxSubmitPU=1 — at most 1 INTR job in queue at a time.
# Use INTR for the first job (amplify, gets immediate scheduling) and
# standard QOS for all others (they queue together without blocking each other).
INTR="scripts/hpc_templates/polychord_intr.sbatch"
STD="scripts/hpc_templates/polychord_standard.sbatch"

echo "=== Stage A: Dark-Photon-Plasma (T1) — paired amplify + suppress ==="

# Stage A amplification: posterior concentrated where P_max is largest
bash scripts/hpc_shuttle.sh submit \
  --template "$INTR" \
  --nodes 1 --ntasks 16 --time 01:00:00 --name "dp_amplify" \
  --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/dp_amplify_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/dark_photon_plasma.json \
  --param kappa=1.0 --param B0=0.01 \
  --prior 'mA2=log_uniform:0.001:1.0' \
  --prior 'deltam=uniform:-0.5:0.5' \
  --prior 'xi=log_uniform:0.05:20.0' \
  --prior 'alpha3=log_uniform:0.001:0.5' \
  --likelihood 'P_max:maximize' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 400 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_amplify.png"

# Stage A suppression: posterior concentrated where P_max is smallest (logL = -P_max)
bash scripts/hpc_shuttle.sh submit \
  --template "$STD" \
  --nodes 1 --ntasks 16 --time 01:00:00 --name "dp_suppress" \
  --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/dp_suppress_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/dark_photon_plasma.json \
  --param kappa=1.0 --param B0=0.01 \
  --prior 'mA2=log_uniform:0.001:1.0' \
  --prior 'deltam=uniform:-0.5:0.5' \
  --prior 'xi=log_uniform:0.05:20.0' \
  --prior 'alpha3=log_uniform:0.001:0.5' \
  --likelihood 'P_max:minimize' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 400 \
  --num-repeats 5 --precision-criterion 0.01 \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_suppress.png"

echo ""
echo "=== Stage B: Einstein-Cartan null (T2) ==="

# Stage B: single amplification run — if null, suppression run is redundant
bash scripts/hpc_shuttle.sh submit \
  --template "$STD" \
  --nodes 1 --ntasks 16 --time 01:00:00 --name "pgt_null" \
  --cmd "OUTPUT_DIR=\${TIDAL_ROOT}/results/pgt_null_\${SLURM_JOB_ID}; mkdir -p \${OUTPUT_DIR}; \
\${MPIRUN_PREFIX} tidal sample examples/data/torsion_gertsenshtein_b5_zero.json \
  --param kappa=1.0 --param B0=0.01 \
  --prior 'alpha1=uniform:-2:2' \
  --prior 'alpha2=uniform:-2:2' \
  --prior 'alpha3=log_uniform:0.01:4' \
  --likelihood 'P_max:maximize' \
  --baseline-formula 'sin(kappa*B0*t_end/2)**2' \
  --method nested --sampler polychord --nlive 300 \
  --num-repeats 4 --precision-criterion 0.01 --no-clustering \
  --grid-shape 64 --bounds 0:50 --periodic \
  --ic plane-wave --ic-component h_5 --ic-wavevector 1.0 --ic-amplitude 1e-2 \
  --source h_5 --target a_1 --snapshots 3 --t-end 10.0 \
  --output \${OUTPUT_DIR} --analyze && \
tidal plot \${OUTPUT_DIR} --type corner --output \${OUTPUT_DIR}/corner_null.png"

echo ""
echo "=== All 3 jobs submitted. ==="
echo "Record job IDs above in CAMPAIGN.md HPC Job Log table."
