#!/usr/bin/env bash
# B.6 — D1 amp v3 bounds-dependence cross-check at L=75 (INTR)
#
# v2 Phase 6.J landed L=75 amp at log Z = +0.842 (Phase 6.G L=100 = +0.679,
# Phase 6.G L=50 = −2.091). Non-monotonic in L; L=75 was the local maximum
# in the v2 architecture. This v3 re-run answers: does the v3 architecture
# preserve the L-dependence pattern, or does the wider posterior support
# wash it out?
#
# Compare against:
#   - v2 reference: 29019705 (Phase 6.J L=75, log Z=+0.842)
#   - v3 L=100 reference: 29149987 (B.1 D1 amp smoke, log Z=+13.29)
#
# Settings: matches B.1 (29149987) except --bounds=0:75 and the wavevector
# adjusted to 2π/75 = 0.0837758... (fundamental mode for the L=75 box).
#
# v3 architecture: probe metadata only; soft floor logL = −100 + Normal(0,1)
# for sim/NaN/exception; arctan_uniform:-89:89 for sign-symmetric; α₃
# log_uniform:1e-3:1e3 (kinetic positivity per v3.1, see GH #358).
#
# See docs/V3_PHASE_TRACKER.md B.6 row.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_amp_L75_v3 --ntasks 32 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 128 --bounds=0:75 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.08377580409572781 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 300 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_amp_L75_v3'
