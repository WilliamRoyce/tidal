#!/usr/bin/env bash
# v3 D1 amp — Ricci-EM, tachyon-permissive, compactified priors
#
# v3 architecture (post-2026-05-08 supervisor pivot):
#   * Probe is metadata only (no gate); --gated would re-enable v2 behavior.
#   * No Hwang-Noh perturbativity gate; no upper P_max cap.
#   * Soft floor logL = -100 + Normal(0, 1.0) for sim failure / NaN / exception.
#   * Compactified priors via arctan_uniform (full real line, Cauchy density at 0).
# See docs/V3_ARCHITECTURE.md and docs/lagrangian_depruning_audit.md.
#
# Reference v2 chain: hpc_results/28896653 (D1 amp v2 hi-res, log Z=+0.720±0.024).
# v3 expectation: posterior shifts measurably; samples that v2 rejected via
# Hwang-Noh (P_max>0.5) or probe (γ_eff>0.15) now contribute. Headline shifts
# from log Z to per-coupling marginal D_KL.

set -euo pipefail

# PolyChord sampling budget — ndim=4:
#   Landscape (this script): nlive=100  (25×ndim), num_repeats=8 (2×ndim)
#   Publication (done):      nlive=200  (50×ndim) via standard queue job 29189748
#                            (logZ=+13.374±0.060, ESS=5694 — supersedes this script)
#
# NOTE: publication results already captured at nlive=1200 (job 29189748).
# This script runs the landscape orientation pass; use for reference or rerun.

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_amp_v3 --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --soft-floor-noise 1.0 \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 512 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 10 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 100 --num-repeats 8 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_amp_v3'
