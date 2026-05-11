#!/usr/bin/env bash
# v3 D2.1 Barker amp RE-RUN — 5D Barker model, undersampling investigation
#
# First v3 run (29209010, nlive=400, num-repeats=5) converged with only
# ESS=474 samples and log Z≈0, MAP at coupling-edge (β1=-19.8, β2=+16.7) —
# classic undersampling signature for 5D PolyChord (80 live points/dim).
# Chain hadn't walked back from prior boundary.
#
# Re-run with nlive=1200 (3x), num-repeats=10 (2x) — 240 live points/dim,
# more slice-sampling steps per iteration.  Expected wall: 30-50 min INTR.
#
# Decision criterion vs 29209010:
# * If log Z still ≈ 0 and ESS scales linearly with nlive → D2.1 Barker
#   is a genuine null under v3 (flat posterior over wide prior support).
# * If log Z deviates substantially (|ΔlogZ| > 1 nat) or chain finds an
#   interior MAP → 29209010 was undersampled; the upgraded settings
#   become the new D2.1 baseline.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d21_barker_sup_v3_rerun --ntasks 76 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
    --param kappa=1.0 --param B0=0.01 \
    --param delta1=0 --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "chi=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1200 \
    --num-repeats 10 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d21_barker_sup_v3_rerun'
