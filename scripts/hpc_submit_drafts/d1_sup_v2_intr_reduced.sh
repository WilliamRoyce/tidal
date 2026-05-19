#!/usr/bin/env bash
# D1 sup INTR reduced — paired counterpart to d1_amp_v2_intr_reduced
#
# Reduced from v2 hi-res: grid=256 (vs 512), nlive=600 (vs 1200), --qos=intr
# Other settings identical to v2 (snapshots=2, k_IC=2π/100, bounds=0:100, etc.)
#
# Submit AFTER amp INTR completes (INTR MaxSubmitPU=1).
# Expected wall: ~30-50 min (sup ~100 clusters at grid=256).

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_sup_v2_intr_reduced --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=uniform:-1:1" --prior "alpha2=uniform:-2:2" \
    --prior "alpha3=log_uniform:0.05:2" --prior "delta1=uniform:-2:2" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 256 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 600 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_sup_v2_intr_reduced'
