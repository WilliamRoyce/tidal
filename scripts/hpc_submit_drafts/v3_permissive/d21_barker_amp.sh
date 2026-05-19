#!/usr/bin/env bash
# v3 D2.1 Barker amp — sub-theory of T5 (6 free params: β₁₋₃, ξ, δ₁, χ)
#
# Parameter partition matches v2 chain 28607124 exactly.
# Pinned: zeta1, zeta2, zeta3 (these were 0 in v2 too).
# Reference v2: 28607124 (log Z = +0.618 ± 0.001).
# v3 29256858: ✅ converged 27 min at nlive=1800 (over-resolved).
#
# PolyChord sampling budget — ndim=6:
#   Landscape (this script):    nlive=150  (25×ndim), num_repeats=12 (2×ndim)
#   Publication (future rerun): nlive=300 (50×ndim), num_repeats=12

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d21_barker_amp_v3 --ntasks 76 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_general_nonminimal.json \
    --param kappa=1.0 --param B0=0.01 \
    --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --grid-shape 64 --bounds=0:50 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 150 \
    --num-repeats 12 --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d21_barker_amp_v3'
