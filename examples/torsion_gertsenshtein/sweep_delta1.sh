#!/usr/bin/env bash
# 1D Delta1 Sweep: Nonminimal R̃[μν]F Coupling Strength
#
# Sweeps delta1 (nonminimal torsion-EM coupling) from -1.0 to 1.0 with fixed
# all-sector-stable torsion mass parameters (alpha1=0, alpha2=-0.6, alpha3=0.5).
#
# Physics:
#   L = (1/κ²) R̃ + α₁I₁ + α₂I₂ + α₃I₃ + δ₁ R̃_{[μν]} F^{μν} − ¼F²
#
#   delta1 = 0: standard Gertsenshtein (no torsion coupling)
#   delta1 > 0: suppression observed at delta1=0.1 (P=0.032 vs P_EM=0.061)
#   delta1 < 0: to be determined (possibly amplification or also suppression)
#
# Mass stability (all sectors negative = Proca-like):
#   Tensor: (1/2)/κ² + (α₁+α₂) = 0.5 + (-0.6) = -0.1  ✓
#   Axial:  (1/24)/κ² - (1/6)(α₁-α₂) = 0.042 - 0.1 = -0.058  ✓
#   Trace:  -(4/3)/κ² + α₃ + (1/3)(2α₁+α₂) = -1.033  ✓
#
# Active channel: h_5 (h× graviton) ↔ a_1 (a_x photon)
#
# Refs:
#   Shapiro (2002), hep-th/0103093 — non-minimal coupling enumeration
#   Barker (2024), arXiv:2406.12826 — PGT particle spectra
#
# Running:
#   cd examples/torsion_gertsenshtein && bash sweep_delta1.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/torsion_gertsenshtein_nonminimal.json"
OUTPUT="${SCRIPT_DIR}/../data/nonminimal_sweep_delta1"

# Fixed physics parameters (all-sector stable)
KAPPA=1.0
B0=0.01
ALPHA1=0
ALPHA2=-0.6
ALPHA3=0.5

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0
K0=2.0

echo "=== Nonminimal R̃[μν]F: Delta1 Sweep ==="
echo "Fixed: kappa=${KAPPA}, B0=${B0}, alpha1=${ALPHA1}, alpha2=${ALPHA2}, alpha3=${ALPHA3}"
echo "Sweep: delta1 = -1.0 .. 1.0 (20 points)"
echo "Gertsenshtein baseline: P = sin^2(kappa * B0 * t_end / 2) = sin^2(0.25) ~ 0.0612"
echo "Output: ${OUTPUT}"
echo ""

uv run tidal sweep "${SPEC}" \
  --sweep "delta1=-1.0:1.0:20" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --param "alpha1=${ALPHA1}" --param "alpha2=${ALPHA2}" --param "alpha3=${ALPHA3}" \
  --resume \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"

# P_max vs delta1 with sin^2(kappa*B0*t_end/2) baseline
uv run tidal plot "${OUTPUT}" --type sweep \
  --metric P_max \
  --title "Nonminimal R̃[μν]F: P vs delta1 (alpha2=${ALPHA2}, alpha3=${ALPHA3})" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2' \
  --output "${OUTPUT}/sweep_delta1.png" --quiet

# C₀ analysis
uv run python "${SCRIPT_DIR}/../torsion_dark_photon/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}"
echo "Plots:"
echo "  ${OUTPUT}/sweep_delta1.png  [P_max vs delta1 with Gertsenshtein baseline]"
echo "  ${OUTPUT}/analysis_C0.png   [C₀ = P/B₀² with baseline]"
