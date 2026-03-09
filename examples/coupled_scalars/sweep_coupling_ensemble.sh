#!/usr/bin/env bash
# Coupled Scalars — Ensemble IC Perturbation Sweep
#
# Demonstrates ensemble repeat runs with IC perturbation for sensitivity
# analysis. Each parameter point is simulated N times with slightly
# different initial conditions (additive Gaussian noise on top of the
# deterministic Gaussian wave packet).
#
# This quantifies how sensitive the conversion probability is to small
# perturbations in the initial state — a standard technique in
# computational physics for assessing IC sensitivity.
#
# Physics:
#   Same coupling sweep as sweep_coupling.sh, but with 10 replicates
#   per point. IC perturbation scale = 5% of ic-amplitude, so each
#   replicate sees a slightly different Gaussian profile.
#
# Statistical output:
#   results.csv     — all individual replicate rows (80 total = 8 x 10)
#   results_agg.csv — one row per parameter point with mean/std/SEM/CI
#   Plots show shaded ±1σ error bands around the mean curve.
#
# References:
#   Palmer, T.N. et al. (1993) "Ensemble prediction", ECMWF Tech Memo 188.
#   Smith, R.C. (2013) Uncertainty Quantification, SIAM. Ch. 3.
#
# Running:
#   cd examples/coupled_scalars && bash sweep_coupling_ensemble.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_ensemble

echo "=== Coupled Scalars — IC Perturbation Ensemble ==="
echo "    gCpl = 0.1 to 1.8 (8 points), 10 replicates per point"
echo "    IC perturbation: 5% of amplitude (Gaussian noise on field slots)"
echo ""

tidal sweep ../data/coupled_scalars.json \
  --sweep "gCpl=0.1:1.8:8" \
  --param mPhi2=1.0 --param mChi2=4.0 \
  --measure peak_conversion,conservation \
  --source phi_0 --target chi_0 \
  --n-replicates 10 --base-seed 42 \
  --ic gaussian --ic-perturbation 0.05 \
  --ic-amplitude 0.1 --ic-width 3.0 --ic-center 20.0 --ic-component phi_0 \
  --grid-shape 128 \
  --bounds 0:160 \
  --periodic \
  --t-end 60.0 \
  --output "$OUT"

echo ""

# Plot with error bands (auto-detected from replicate data)
tidal plot "$OUT" --type sweep \
  --metric P_max \
  --title "Coupling sweep with IC perturbation (N=10)" \
  --output "$OUT/sweep.png" --quiet

echo ""
echo "=== Ensemble sweep complete ==="
echo "Results: $OUT/"
echo "  results.csv      — all 80 replicate rows"
echo "  results_agg.csv  — 8 aggregated rows (mean/std/SEM/CI)"
echo "  sweep.png        — P_max vs gCpl with ±1σ error bands"
