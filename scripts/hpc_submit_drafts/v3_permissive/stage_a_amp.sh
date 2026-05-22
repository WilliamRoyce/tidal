#!/usr/bin/env bash
# v3 Stage A amp — Dark-Photon-Plasma, tachyon-permissive, compactified priors
#
# v3 architecture: see docs/V3_ARCHITECTURE.md.
# Stage A v2 result was a clean null (log Z = -0.073 ± 0.007, D_KL=0.024 nats).
# v3 re-run: do tachyon-permissive sampling and wider priors change the
# verdict, or is the null robust?  α₃ MAP at v2 ≈ 0.001 was far from any
# probe boundary, so we expect minimal shift — but the test is the point.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name stage_a_amp_v3 --ntasks 76 \
  --cmd 'tidal sample examples/data/dark_photon_plasma.json \
    --prior "mA2=log_uniform:1e-3:1e3" \
    --prior "deltam=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/stage_a_amp_v3'
