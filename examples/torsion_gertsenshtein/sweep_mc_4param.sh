#!/usr/bin/env bash
# 4-Parameter Monte Carlo Sweep: Nonminimal R̃[μν]F Torsion Model
#
# Latin Hypercube sampling over (α₁, α₂, α₃, δ₁) with fixed κ=1, B₀=0.001.
# Maps amplification landscape of the Gertsenshtein effect in PGT with
# non-minimal torsion-EM coupling.
#
# Physics:
#   L = (1/κ²) R̃ + α₁I₁ + α₂I₂ + α₃I₃ + δ₁ R̃_{[μν]} F^{μν} − ¼F²
#
# Stability constraints (all sectors Proca-like):
#   Tensor: α₁ + α₂ < -0.5
#   Axial:  α₁ - α₂ > 0.25
#   Trace:  -(4/3) + α₃ + (1/3)(2α₁ + α₂) < 0
#
# Unstable points will diverge (caught by divergence guard); stable points
# give P_final from which we compute:
#   C₀ = P_final / B₀²           (conversion coefficient)
#   A  = C₀ / C_EM               (amplification factor vs standard GR)
#
# Active channel: h_5 (h× graviton) ↔ a_1 (a_x photon)
#
# Running:
#   cd examples/torsion_gertsenshtein && bash sweep_mc_4param.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/torsion_gertsenshtein_nonminimal.json"
OUTPUT="${SCRIPT_DIR}/../data/nonminimal_mc_4param"

# Fixed physics parameters
KAPPA=1.0
B0=0.001  # Small B₀ for linear regime (P << 1)

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0
K0=2.0
N_SAMPLES=500

echo "=== Nonminimal R̃[μν]F: 4-Parameter LHS Sweep ==="
echo "Fixed: kappa=${KAPPA}, B0=${B0}"
echo "Swept: alpha1∈[-2,2], alpha2∈[-3,0], alpha3∈[0,5], delta1∈[-1,1]"
echo "Strategy: Latin Hypercube, ${N_SAMPLES} samples"
echo "Output: ${OUTPUT}"
echo ""

uv run tidal sweep "${SPEC}" \
  --sweep "alpha1=-2:2:10" --sweep "alpha2=-3:0:10" \
  --sweep "alpha3=0:5:10" --sweep "delta1=-1:1:10" \
  --sweep-strategy latin_hypercube --n-samples "${N_SAMPLES}" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --resume --force \
  --output "${OUTPUT}"

echo ""
echo "--- Generating analysis plots ---"

# Tornado chart: parameter sensitivity
uv run tidal plot "${OUTPUT}" --type sweep-tornado \
  --metric P_final \
  --output "${OUTPUT}/tornado.png" --quiet

# Parallel coordinates: parameter-metric relationships
uv run tidal plot "${OUTPUT}" --type sweep-parallel \
  --metric P_final \
  --output "${OUTPUT}/parallel.png" --quiet

# Pairwise scatter: interaction structure
uv run tidal plot "${OUTPUT}" --type sweep-scatter \
  --metric P_final \
  --output "${OUTPUT}/scatter.png" --quiet

# C₀ amplification analysis
uv run python "${SCRIPT_DIR}/../torsion_dark_photon/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}/results.csv"
echo "Plots:"
echo "  ${OUTPUT}/tornado.png     [Parameter sensitivity]"
echo "  ${OUTPUT}/parallel.png    [Parallel coordinates]"
echo "  ${OUTPUT}/scatter.png     [Pairwise scatter matrix]"
echo "  ${OUTPUT}/analysis_C0.png [C₀ amplification analysis]"
