#!/usr/bin/env bash
# v3 D2.0 Bahamonde amp — sub-theory of T5 (5 free params: β₁₋₃, ξ, δ₁)
#
# v3 architecture: see docs/V3_ARCHITECTURE.md.
# Parameter partition matches v2 chain 28598736 exactly (verified
# against actual inference.json, not just the legacy submission script).
# Pinned: chi, zeta1, zeta2, zeta3 (these were 0 in v2 too).
#
# Reference v2: 28598736 (log Z = +0.616 ± 0.001, ESS=877).
# v2 priors were tight (xi=log[0.01,10], delta1=uniform[-0.025,0.025]);
# v3 widens them via the architecture-standard compactification.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d20_bahamonde_amp_v3 --ntasks 76 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
    --param kappa=1.0 --param B0=0.01 \
    --param chi=0 --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 1500 \
    --num-repeats 5 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d20_bahamonde_amp_v3'
