#!/usr/bin/env bash
# Phase 6.M — D1 amp at t_end=3 (perturbativity-window cross-check)
#
# Phase 6.L found that all v2 amp MAPs (τ=0.15) sit at the boundary
# attractor γ_eff ≈ 0.13–0.15, so the chain's A_chain(t_end=10) factorises as
#
#     A_chain(t_end) = A_static × exp(γ_eff · t_end)
#                       │            │
#                       │            └ residual tachyonic growth (contamination)
#                       └ light-mediator (publishable)
#
# At γ_eff = 0.15, t_end=10 → exp(γ·t_end) ≈ 4.5× residual factor. Shortening
# to t_end=3 brings this to exp(0.45) ≈ 1.6× — A_chain is then ≈A_static within
# 60% and the perturbativity ratio A(2t)/A(t) ≈ 1.5 sits cleanly inside the
# [0.5, 2.0] pass band even at the τ=0.15 probe boundary.
#
# Probe still runs at t_test=20 (long-horizon stability screen). Asymmetry
# probe>>chain is by design — probe catches slow tachyons; simulation reads
# physics in a clean linear window.
#
# Expected results (predicted from Phase 6.L analysis):
#   * log Z lower than v2 t_end=10 (A_static ≈ 10, baseline ~11× smaller, but
#     the baseline cancels in the metric ratio so log Z penalty comes only
#     from the smaller exp(γ·t) factor: Δlog Z ≈ −log(exp(γ·7)) ≈ −1 nat)
#   * MAP roughly at same (α₁,α₂,α₃,δ₁) (boundary attractor still exists,
#     just less amplified)
#   * A_chain at MAP ≈ 15 (vs 38 at t_end=10)
#   * A_static = A_chain / exp(γ_eff · 3) ≈ 10 — the headline number
#
# Settings: matched to v2 INTR-reduced (28982029, grid=256/nlive=600) for
# direct comparison. Only t_end and snapshots count are changed.
#
# Reference: examples/data/d1_perturbativity_check_v2/SUMMARY.md
# Pairs with: 28982029 (v2 amp INTR-reduced at t_end=10)

set -euo pipefail

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name d1_amp_tend3_intr --ntasks 76 \
  --cmd 'tidal sample examples/data/torsion_gertsenshtein_nonminimal.json \
    --prior "alpha1=uniform:-1:1" --prior "alpha2=uniform:-2:2" \
    --prior "alpha3=log_uniform:0.05:2" --prior "delta1=uniform:-2:2" \
    --likelihood "P_max:maximize" \
    --baseline-formula "sin(kappa*B0*t_end/2)**2" \
    --param kappa=1.0 --param B0=0.01 \
    --grid-shape 256 --bounds=0:100 --periodic \
    --ic plane-wave --ic-component h_5 --ic-wavevector 0.06283185307179587 --ic-amplitude 1e-2 \
    --t-end 3 --snapshots 2 \
    --measure conversion,peak_conversion --source h_5 --target a_1 \
    --method nested --sampler polychord --nlive 600 \
    --precision-criterion 0.01 \
    --output ${RESULTS_DIR}/d1_amp_tend3_intr'
