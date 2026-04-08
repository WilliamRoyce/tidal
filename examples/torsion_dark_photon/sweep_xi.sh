#!/usr/bin/env bash
# 1D Xi Sweep: Torsion Amplification vs Kinetic Coupling Strength
#
# Sweeps xi (torsion kinetic term coefficient) from 0.01 to 0.5 with fixed
# alpha=0.5 and deltam=0.1. Measures conversion probability P_max at each xi.
#
# Physics:
#   L = (1/κ²)R̃ + α I₃ - ¼ξ Ftorsion² + δ F·Ftorsion - ¼F²
#
#   xi → 0: torsion kinetic term vanishes → torsion becomes algebraically
#     constrained (time_order=0) → decouples from dynamics → exact Gertsenshtein.
#   xi > 0: torsion propagates as dark photon, kinetic mixing δ F·Ftorsion
#     introduces three-field graviton↔photon↔torsion mixing → amplification.
#
# Expected output:
#   P_max(xi=0⁺) ≈ P_Gertsenshtein = sin²(κ B0 t_end / 2)
#               = sin²(1.0 × 0.01 × 50 / 2) = sin²(0.25) ≈ 0.0617
#
#   P_max(xi=0.1) ≈ 7.24e-4 * B0² * (t_end/D)²  [torsion-enhanced]
#   ... actually P_max rises above baseline for xi > 0 when alpha, deltam > 0.
#
# The overlay `sin(kappa * B0 * t_end / 2)**2` is the xi=0 Gertsenshtein
# baseline. Points above it show amplification from torsion.
#
# Note: At small B0=0.01 we are in the small-P regime where C0 = P/B0^2
# is parameter-independent of B0. The analytical baseline is P ≈ 0.0617.
#
# Active channel: h_5 (h× graviton) ↔ a_1 (a_x photon)
#   coupling: -B0*kappa^2 * gradient_z in h_5 EOM
#
# Refs:
#   Barker et al. (2024), arXiv:2406.12826 — eWGT / propagating vector torsion
#   Blagojević & Cvetković (2019), arXiv:1812.02675 — ghost/tachyon-free windows
#
# Running:
#   cd examples/torsion_dark_photon && bash sweep_xi.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/torsion_dark_photon.json"
OUTPUT="${SCRIPT_DIR}/../data/torsion_dark_photon_sweep_xi"

# Fixed physics parameters
KAPPA=1.0
B0=0.01        # small-P regime
K0=2.0         # plane-wave wavevector
ALPHA=0.5      # torsion mass (trace sector)
DELTAM=0.1     # EM-torsion kinetic mixing

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

echo "=== Dark Photon Torsion: Xi Sweep ==="
echo "Fixed: kappa=${KAPPA}, B0=${B0}, alpha=${ALPHA}, deltam=${DELTAM}"
echo "Sweep: xi = 0.01 .. 0.5 (20 points)"
echo "Gertsenshtein baseline: P = sin^2(kappa * B0 * t_end / 2) = sin^2(0.25) ~ 0.0617"
echo "Output: ${OUTPUT}"
echo ""

# Note: --measure peak_conversion only (no conservation).
# The total Hamiltonian is NOT a valid conservation diagnostic for this constrained
# gauge theory (Dirac-Bergmann). Only the h_5↔a_1 physical sector has a
# well-defined conserved energy. Use P_max convergence and run_status instead.
uv run tidal sweep "${SPEC}" \
  --sweep "xi=0.01:0.5:20" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --param "alpha=${ALPHA}" --param "deltam=${DELTAM}" \
  --no-require-stable \
  --resume \
  --output "${OUTPUT}"
# --no-require-stable: the pre-simulation mass check flags R̃ tensor/axial
# eigenvalues as tachyonic, but these modes have ZERO physical coupling
# (100% trace-aligned). The modal solver's _suppress_tachyonic_noise() (#222)
# freezes these modes during evolution, preventing numerical noise growth.

echo ""
echo "--- Generating plots ---"

# P_max vs xi with sin^2(kappa*B0*t_end/2) as the xi=0 Gertsenshtein baseline.
# Points above the dashed line = torsion amplification.
uv run tidal plot "${OUTPUT}" --type sweep \
  --metric P_max \
  --title "Torsion amplification vs xi (alpha=${ALPHA}, deltam=${DELTAM}, B0=${B0})" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2' \
  --output "${OUTPUT}/sweep_xi.png" --quiet

# C₀ = P/B₀² analysis (primary metric per project_sweep_measurement_strategy.md)
uv run python "${SCRIPT_DIR}/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}"
echo "Plots:"
echo "  ${OUTPUT}/sweep_xi.png     [raw P_max vs xi]"
echo "  ${OUTPUT}/analysis_C0.png  [C₀ = P/B₀² with Gertsenshtein baseline]"
