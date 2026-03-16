#!/usr/bin/env bash
# Gertsenshtein Effective — Parameter Noise Uncertainty Quantification
#
# Demonstrates uncertainty propagation via parameter noise. The plasma
# frequency omegaP2 is treated as uncertain: at each replicate, its value
# is drawn from N(omegaP2_nominal, σ²) where σ = 0.01.
#
# This quantifies how sensitive the conversion probability is to
# uncertainties in the plasma environment — relevant for astrophysical
# applications where the plasma density is not precisely known.
#
# Ref: Smith, R.C. (2013) Uncertainty Quantification, SIAM. Ch. 3-5.
#      Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237.
#
# Running:
#   cd examples/coupled_scalars && bash sweep_param_uq.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_param_uq

echo "=== Gertsenshtein Effective — Parameter Noise UQ ==="
echo "    B0 = 0.01 to 0.25 (8 points), 5 replicates per point"
echo "    Parameter noise: omegaP2 ~ N(0.0, 0.01²)"
echo ""

tidal sweep ../data/coupled_scalars.json \
  --sweep "B0=0.01:0.25:8" \
  --param kappa=1.0 --param omegaP2=0.0 --param mg2=0.0 \
  --measure conversion,conservation \
  --source h_0 --target a_0 \
  --scheme cvode \
  --n-replicates 5 --base-seed 42 \
  --param-noise "omegaP2=0.01" \
  --ic plane-wave --ic-amplitude 0.1 --ic-wavevector 2.0 --ic-component h_0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --t-end 50.0 \
  --output "$OUT"

echo ""

tidal plot "$OUT" --type sweep \
  --metric P_final \
  --title "B0 sweep with plasma uncertainty (σ_omegaP2 = 0.01)" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2' \
  --output "$OUT/sweep.png" --quiet

echo ""
echo "=== UQ sweep complete ==="
echo "Results: $OUT/"
echo "  sweep.png — P_final vs B0 with ±1σ uncertainty bands"
