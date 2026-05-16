#!/usr/bin/env bash
# v3 D2.3 full T5 sup — all 9 nonminimal couplings free
#
# Reference v2: 28742035 (log Z = -0.615 ± 0.001, std 6h, 1:00 wall).
# max_ndead=2000: very conservative for 9p sup (most complex chain); ~15 min + post-processing.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d23_full_sup_v3 --ntasks 76 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
    --param kappa=1.0 --param B0=0.01 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1800 \
    --num-repeats 5 --precision-criterion 0.01 \
    --max-ndead 2000 \
    --output ${RESULTS_DIR}/d23_full_sup_v3'
