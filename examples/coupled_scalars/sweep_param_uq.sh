#!/usr/bin/env bash
# Coupled Scalars — Parameter Noise Uncertainty Quantification
#
# Demonstrates uncertainty propagation via parameter noise. The mass
# parameter m2 is treated as uncertain: at each replicate, its value
# is drawn from N(m2_nominal, σ²) where σ = 0.1 (5% of m2=2.0).
#
# This is the standard forward uncertainty propagation approach:
# given uncertain inputs, characterize the output distribution.
#
# Physics:
#   Coupling sweep with uncertain mass ratio. The mass gap Δm² = m1² - m2²
#   controls the Rabi oscillation frequency, so mass uncertainty directly
#   affects the conversion probability.
#
# Statistical output:
#   results.csv     — all replicate rows with realized m2 values
#   results_agg.csv — aggregated rows with {metric}_mean, _std, _sem, _ci95
#   Columns include m2_nominal (grid value) and m2 (realized value with noise)
#
# References:
#   Smith, R.C. (2013) Uncertainty Quantification, SIAM. Ch. 3-5.
#   Saltelli, A. et al. (2008) Global Sensitivity Analysis, Wiley. Ch. 1.
#
# Running:
#   cd examples/coupled_scalars && bash sweep_param_uq.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_param_uq

echo "=== Coupled Scalars — Parameter Noise UQ ==="
echo "    gCpl = 0.1 to 1.8 (8 points), 10 replicates per point"
echo "    Parameter noise: mChi2 ~ N(4.0, 0.1²)"
echo ""

tidal sweep ../data/coupled_scalars.json \
  --sweep "gCpl=0.1:1.8:8" \
  --param mPhi2=1.0 --param mChi2=4.0 \
  --measure peak_conversion,conservation \
  --source phi_0 --target chi_0 \
  --n-replicates 10 --base-seed 42 \
  --param-noise "mChi2=0.1" \
  --ic gaussian --ic-amplitude 0.1 --ic-width 3.0 --ic-center 20.0 --ic-component phi_0 \
  --grid-shape 128 \
  --bounds 0:160 \
  --periodic \
  --t-end 60.0 \
  --output "$OUT"

echo ""

# Plot with error bands
tidal plot "$OUT" --type sweep \
  --metric P_max \
  --title "Coupling sweep with mass uncertainty (σ_mChi2 = 0.1)" \
  --output "$OUT/sweep.png" --quiet

echo ""
echo "=== UQ sweep complete ==="
echo "Results: $OUT/"
echo "  results.csv      — all 80 replicate rows (with mChi2_nominal column)"
echo "  results_agg.csv  — 8 aggregated rows (mean/std/SEM/CI)"
echo "  sweep.png        — P_max vs gCpl with ±1σ uncertainty bands"
