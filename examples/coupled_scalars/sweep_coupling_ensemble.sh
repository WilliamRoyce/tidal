#!/usr/bin/env bash
# Gertsenshtein Effective — Ensemble IC Perturbation Sweep
#
# Demonstrates ensemble repeat runs with IC perturbation for uncertainty
# quantification. Each B0 point is simulated N times with slightly different
# initial conditions (additive Gaussian noise on the plane-wave IC).
#
# Uses CVODE scheme (not modal) because modal solver's expm_multiply
# path has genuinely growing modes for certain fields when ICs are perturbed
# away from the physical subspace.
#
# Statistical output:
#   results.csv     — all individual replicate rows
#   results_agg.csv — one row per parameter point with mean/std/SEM/CI
#   Plots show shaded ±1σ error bands around the mean curve.
#
# Ref: Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#
# Running:
#   cd examples/coupled_scalars && bash sweep_coupling_ensemble.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_ensemble

echo "=== Gertsenshtein Effective — IC Perturbation Ensemble ==="
echo "    B0 = 0.01 to 0.25 (8 points), 5 replicates per point"
echo "    IC perturbation: 2% of amplitude"
echo ""

tidal sweep ../data/coupled_scalars.json \
  --sweep "B0=0.01:0.25:8" \
  --param kappa=1.0 --param omegaP2=0.0 --param mg2=0.0 \
  --measure conversion,conservation \
  --source h_0 --target a_0 \
  --scheme cvode \
  --n-replicates 5 --base-seed 42 \
  --ic plane-wave --ic-perturbation 0.02 \
  --ic-amplitude 0.1 --ic-wavevector 2.0 --ic-component h_0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --t-end 50.0 \
  --output "$OUT"

echo ""

tidal plot "$OUT" --type sweep \
  --metric P_final \
  --title "B0 sweep with IC perturbation (N=5)" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2' \
  --output "$OUT/sweep.png" --quiet

echo ""
echo "=== Ensemble sweep complete ==="
echo "Results: $OUT/"
echo "  sweep.png — P_final vs B0 with ±1σ error bands"
