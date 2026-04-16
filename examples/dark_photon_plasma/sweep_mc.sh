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
#   mA2    ∈ [0.001, 500]  — photon effective mass (plasma frequency)
#   mT2    ∈ [0.001, 500]  — dark-photon Proca mass
#   deltam ∈ [-50, 50]     — kinetic mixing (both signs, full range)
#
# Focused bounds for the physically interesting regime:
#   mA² ∈ [0.01, 5] — covers plasma frequencies relevant for
#     astrophysical environments and the resonance condition mA² ≈ mT².
#   mT² ∈ [0.01, 5] — matched to mA² so resonance falls within range.
#   δ ∈ [-2, 2] — covers weak mixing (|δ| ≪ 1), the Holdom critical
#     point (|δ| = √ξ/2 = 0.5 for ξ=1), and the strong-mixing regime
#     (|δ| > 1 where δ² > ξ → kinetic matrix becomes ghost-like).
# Arctan mapping: ~50% of samples in [0, 1], ~50% in [1, 5].
MA2_BOUNDS="0.01:5.0:2"
MT2_BOUNDS="0.01:5.0:2"
DELTA_BOUNDS="-2.0:2.0:2"

echo "=== Plasma Dark Photon: Monte Carlo Sweep ==="
echo "Samples: ${N_SAMPLES}"
echo "Fixed: kappa=${KAPPA}, B0=${B0}, xi=${XI}, k0=${K0}"
echo "Sweep: mA2 ∈ [0.001, 500], mT2 ∈ [0.001, 500], deltam ∈ [-50, 50]"
echo "Output: ${OUTPUT}"
echo ""

tidal sweep "${SPEC}" \
  --sweep "mA2=${MA2_BOUNDS}" \
  --sweep "mT2=${MT2_BOUNDS}" \
  --sweep "deltam=${DELTA_BOUNDS}" \
  --sweep-strategy latin_hypercube \
  --n-samples "${N_SAMPLES}" \
  --measure conversion,peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" --param "xi=${XI}" \
  --parallel "${TIDAL_PARALLEL:-4}" --resume \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"

# Sensitivity: which of (mA2, mT2, deltam) drives conversion?
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
echo "  Resonance condition: mA² ≈ mT² (matched photon/dark-photon dispersions)"
echo "  Holdom triviality: E[t,a] = -2·deltam·mA²/(4·deltam²-xi)"
echo "  Expected P(h→t) ∝ deltam² · mA² at weak mixing"
