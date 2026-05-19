#!/usr/bin/env bash
# D1 amp v2 — publication-quality paired run at tightened τ=0.15 probe
#
# v2 changes vs 28789437/28789439 baseline:
#   (a) probe threshold 0.30 → 0.15 via HEAD commit 98c87d7
#   (b) nlive 400 → 1200
#   (c) grid 256 → 512
#   (d) snapshots default(21) → 2 (user policy)
# All other settings identical to Phase 6.C.2 baseline.
#
# Reference: hpc_results/28789437/ (canonical D1 amp clean, std icelake 6h, 19:24 wall)
# Probe profile: unit-ic-all-k-0.15 (committed pre-Phase-6.G in 98c87d7)
# IC k = 2π/100 = 0.06283185307179587 (fundamental mode of L=100 box)

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_standard.sbatch \
  --name d1_amp_track2_hires_v2 --ntasks 76 --time 08:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=uniform:-1:1" --prior "alpha2=uniform:-2:2" \
    --prior "alpha3=log_uniform:0.05:2" --prior "delta1=uniform:-2:2" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1200 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_amp_track2_hires_v2'
