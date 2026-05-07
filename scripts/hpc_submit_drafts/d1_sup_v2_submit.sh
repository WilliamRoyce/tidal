#!/usr/bin/env bash
# D1 sup v2 — publication-quality paired run at tightened τ=0.15 probe
#
# v2 changes vs 28789439 baseline (Phase 6.C D1 sup clean, std icelake 6h, 2:32 wall, 121 clusters):
#   (a) probe threshold 0.30 → 0.15 via HEAD commit 98c87d7
#   (b) nlive 400 → 1200
#   (c) grid 256 → 512
#   (d) snapshots default(21) → 2 (user policy)
# All other settings identical to Phase 6.C baseline.
#
# Reference: hpc_results/28789439/ (canonical D1 sup clean, std icelake 6h)
# Pairs with d1_amp_track2_hires_v2 — identical params except --likelihood mode.
#
# Cost estimate vs 28789439 (2:32 wall):
#   snapshots 21→2: ~10× faster per call
#   grid 256→512: ~2× slower per call
#   nlive 400→1200 + cluster growth: ~3-4× more dead points
#   τ=0.15 probe accepts more borderline samples → more multi-modal sup surface
#   Empirical: intr_reduced (nlive=600/grid=256) timed out at 35 clusters in 1h (28983285)
#   At v2 settings ~4-6× more expensive per eval → 12h wall (job 28982018 at 6h cancelled)

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_standard.sbatch \
  --name d1_sup_track2_hires_v2 --ntasks 76 --time 12:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=uniform:-1:1" --prior "alpha2=uniform:-2:2" \
    --prior "alpha3=log_uniform:0.05:2" --prior "delta1=uniform:-2:2" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1200 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_sup_track2_hires_v2'
