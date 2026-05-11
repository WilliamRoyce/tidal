#!/usr/bin/env bash
# v3 D1 amp B.0a JOINT-PRIOR SMOKE — head-to-head vs B.1 (29149987)
#
# Same INTR settings (grid=128/nlive=300/t_end=10) as B.1 but uses the
# cubed-sphere joint prior for the 3 sign-symmetric couplings (α₁, α₂, δ₁).
# α₃ stays on per-param log_uniform (positive-definite kinetic, doesn't
# belong on a sign-symmetric sphere).
#
# Joint prior:
#   names=alpha1,alpha2,delta1 → vector c ∈ R^3 = r * θ̂ with θ̂ ∈ S^2
#   M=1: no sub-tile subdivisions (full face)
#   face=1: face 1+ (positive x_1 hemisphere)
#   sub=1_1: only sub-tile when M=1
#   r_lo=1e-3, r_hi=1e3: log_uniform magnitude scale (matches α₃ scale)
#
# Decision criterion (post-run): if ESS ≥ 70% of B.1's 1410 AND walltime
# ≤ 2× B.1's 23 min AND MAP agrees with B.1 within 1σ → adopt joint prior
# for B.5 D2.0–D2.3 (where higher-N spheres matter most).
#
# Output: ${RESULTS_DIR}/d1_amp_v3_jointprior_smoke/01p_tile1_1/ (per-tile
# directory convention).  Pull with --src to the tile sub-directory.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_amp_v3_joint_smoke --ntasks 32 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --joint-prior "names=alpha1,alpha2,delta1;type=cubed_sphere;M=1;face=1;sub=1_1;r_lo=1e-3;r_hi=1e3" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 128 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 300 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_amp_v3_joint_smoke'
