#!/usr/bin/env bash
# v3 D1 amp SMOKE (Phase B.1) — mid-res INTR, tachyon-permissive, compactified priors
#
# Phase B.1 verification criteria (see docs/V3_PHASE_TRACKER.md):
#   * Chain converges within 1 h INTR slot.
#   * Posterior shape includes samples that v2 hard-rejected (probe gate / Hwang-Noh).
#   * Compare MAP to v2 MAP (hpc_results/28982029, 28896653); expect measurable shift.
#   * Four run_status fractions logged; any > 30% → investigate before Phase B.2.
#
# v3 architecture (post-2026-05-08 supervisor pivot):
#   * Probe is metadata only (no gate); --gated would re-enable v2 behaviour.
#   * No Hwang-Noh perturbativity gate; no upper P_max cap.
#   * Soft floor logL = -100 + Normal(0, 1.0) for sim failure / NaN / exception.
#   * Compactified priors via arctan_uniform (full real line, Cauchy density at 0).
# See docs/V3_ARCHITECTURE.md.

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_amp_v3_smoke --ntasks 32 --time 01:00:00 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=arctan_uniform:-89:89" \
    --prior "alpha2=arctan_uniform:-89:89" \
    --prior "alpha3=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
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
    --output ${RESULTS_DIR}/d1_amp_v3_smoke'
