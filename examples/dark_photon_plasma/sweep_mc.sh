#!/usr/bin/env bash
# Monte Carlo Parameter Sweep: Plasma Dark Photon Model
#
# Sweeps over 3 parameters (mA2, mT2, deltam) using Latin Hypercube
# sampling with arctan parameterization to cover a wide range from
# near-zero to hundreds, with dense sampling near the origin.
#
# Physics:
#   L = (1/κ²) R − ¼F² − ½mA²·a² − ¼ξFt² − ½mT²·t² + δ·F·Ft
#
# The photon effective mass mA² (plasma frequency proxy) breaks the
# Holdom triviality of pure kinetic mixing: in vacuum (mA²=0), the
# photon is an exact eigenmode and no dark-photon conversion occurs.
# In plasma (mA²>0), the effective coupling E[t,a] = -2δ·mA²/(4δ²-ξ)
# is non-zero, enabling genuine h → γ(plasma) → γ' conversion.
#
# The Raffelt-Stodolsky formula gives:
#   P(γ→γ') = ε²·sin²((mT²−mA²)·t/(4ω))
# with resonance at mT² ≈ mA² (matched dispersions).
#
# Metrics:
#   A_total = P(h→a+t) / P_GR     — total amplification vs Gertsenshtein
#   A_dark  = P_full / P_plasma    — dark-photon amplification vs plasma-only
#
# Refs:
#   Holdom (1986), Phys. Lett. B 166, 196 — kinetic mixing triviality
#   An, Pospelov, Pradler (2013), arXiv:1302.3884 — conversion formulas
#   Raffelt & Stodolsky (1988), Phys. Rev. D 37, 1237 — two-state mixing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/dark_photon_plasma.json"
OUTPUT="${SCRIPT_DIR}/../data/dark_photon_plasma_mc"

# Number of Monte Carlo samples
N_SAMPLES="${1:-1000}"

# Fixed physics parameters
KAPPA=1.0
B0=0.01        # small-P regime
K0=2.0         # GW frequency (plane-wave wavevector)
XI=1.0         # canonical kinetic normalisation

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

# Sweep parameter bounds (arctan-mapped for wide coverage):
#   mA2    ∈ [0.01, 5]    — photon effective mass (plasma frequency)
#   alpha3 ∈ [0.005, 2.5] — PGT I3 coefficient; mass^2 = 2*alpha3 ∈ [0.01, 5]
#                            (equiv. to old mT2 ∈ [0.01, 5] in FV model)
#   deltam ∈ [-2, 2]      — kinetic mixing (both signs)
#   xi     ∈ [0.1, 5]     — torsion kinetic strength
#
# Resonance condition: 2*alpha3 ≈ mA2 (dark-photon mass² ≈ photon plasma mass²).
# Ghost boundary: |δ| < √ξ/2 (kinetic matrix K_kin = [[1,-2δ],[-2δ,ξ]] positive-definite).
# Arctan mapping: ~50% of samples in [0, 1], ~50% in [1, 5].
MA2_BOUNDS="0.01:5.0:2"
ALPHA3_BOUNDS="0.005:2.5:2"
DELTA_BOUNDS="-2.0:2.0:2"
XI_BOUNDS="0.1:5.0:2"

echo "=== Plasma Dark Photon: 4D Monte Carlo Sweep ==="
echo "Samples: ${N_SAMPLES}"
echo "Fixed: kappa=${KAPPA}, B0=${B0}, k0=${K0}"
echo "Sweep: mA2 ∈ [0.01, 5], alpha3 ∈ [0.005, 2.5] (mass^2=2*alpha3), deltam ∈ [-2, 2], xi ∈ [0.1, 5]"
echo "Ghost boundary: |deltam| < sqrt(xi)/2"
echo "Output: ${OUTPUT}"
echo ""

tidal sweep "${SPEC}" \
  --sweep "mA2=${MA2_BOUNDS}" \
  --sweep "alpha3=${ALPHA3_BOUNDS}" \
  --sweep "deltam=${DELTA_BOUNDS}" \
  --sweep "xi=${XI_BOUNDS}" \
  --sweep-strategy latin_hypercube \
  --n-samples "${N_SAMPLES}" \
  --measure conversion,peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --parallel "${TIDAL_PARALLEL:-4}" --resume \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"

# Sensitivity: which of (mA2, mT2, deltam, xi) drives conversion?
tidal plot "${OUTPUT}" --type sweep-tornado \
  --metric P_max \
  --title "Plasma dark photon: parameter sensitivity (P_max)" \
  --output "${OUTPUT}/plot_tornado.png" --quiet

# Scatter matrix: all parameter pairs coloured by P_max
tidal plot "${OUTPUT}" --type sweep-scatter \
  --metric P_max \
  --title "Plasma dark photon: parameter space (coloured by P_max)" \
  --output "${OUTPUT}/plot_scatter.png" --quiet

# Parallel coordinates: multi-parameter trends
tidal plot "${OUTPUT}" --type sweep-parallel \
  --metric P_max \
  --title "Plasma dark photon: parallel coordinates (P_max)" \
  --output "${OUTPUT}/plot_parallel.png" --quiet

# Analysis with stability filter + paired baseline
python "${SCRIPT_DIR}/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}"
echo "Plots:"
echo "  ${OUTPUT}/plot_tornado.png   [parameter sensitivity ranking]"
echo "  ${OUTPUT}/plot_scatter.png   [pairwise scatter, coloured by P_max]"
echo "  ${OUTPUT}/plot_parallel.png  [parallel coordinates]"
echo ""
echo "Physics notes:"
echo "  Resonance condition: 2*alpha3 ≈ mA2 (dark-photon mass² = photon plasma mass²)"
echo "  Holdom triviality: E[t,a] = -2·deltam·mA²/(4·deltam²-xi)"
echo "  Expected P(h→a) ∝ deltam² · mA² at weak mixing"
