#!/usr/bin/env bash
# v3 D1 sup — Ricci-EM, tachyon-permissive, compactified priors
#
# v3 architecture: see docs/V3_ARCHITECTURE.md and the d1_amp.sh header.
# Sup chains: P_max:minimize captures integrated suppression (IC overlap, phase
# mismatch, destructive interference), which the eigenvalue-based γ_conversion
# cannot.  Sim runs every sample.
#
# Reference v2 chain: hpc_results/28519675 (D1 sup v5, log Z=+15.92±0.13,
# joint D_KL=8.91 nats; deepest suppression valley reaches A ≈ 4×10⁻¹²).
# v3 expectation: chain explores tachyonic regions where v2 would have
# rejected — new sup-relevant structure may emerge ("where in parameter
# space tachyons FAIL to develop while conversion stays low").

set -euo pipefail

# PolyChord sampling budget — ndim=4:
#   Landscape (this script): nlive=100  (25×ndim), num_repeats=8 (2×ndim)
#   Publication (done):      nlive=200  (50×ndim) via standard queue job 29189761
#                            (logZ=+11.395±0.097, ESS=5224 — supersedes this script)

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_sup_v3 --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_sup_v3'
