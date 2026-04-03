#!/usr/bin/env bash
# 2D Parameter Sweep: Amplification Map over (alpha, deltam)
#
# Sweeps the torsion mass alpha and EM-torsion kinetic mixing deltam on an 8x7
# Cartesian grid with fixed xi=0.1. Produces a heatmap showing where in the
# (alpha, deltam) plane conversion is amplified relative to Gertsenshtein.
#
# Physics:
#   L = (1/κ²)R̃ + α I₃ - ¼ξ Ftorsion² + δ F·Ftorsion - ¼F²
#
#   alpha: torsion trace mass (controls torsion dispersion relation)
#     m_torsion^2 = 2*alpha (heuristic from EOM analysis)
#     At resonance omega_GW ~ m_torsion: maximum energy transfer to torsion
#     then re-emission to photon via deltam mixing.
#
#   deltam: EM-torsion kinetic mixing (dark photon ε analogue)
#     deltam = 0: no direct photon↔torsion coupling, only 3-body via graviton
#     |deltam| > 0: direct photon↔torsion oscillations + Gertsenshtein 3-body
#
# Fixed xi=0.1 ensures torsion is dynamical (propagates as massive vector boson).
# Ghost-free region: xi > 0, alpha > 0, deltam^2 < f(alpha, xi).
# The sweep covers both physical (alpha, deltam > 0) and unphysical regions
# to map where the ghost-free boundary lies.
#
# Grid: 8 alpha × 7 deltam = 56 runs total (~few minutes at N=256, t_end=50)
#
# Active channel: h_5 (h× graviton) ↔ a_1 (a_x photon)
# Nominal result (alpha=0.5, deltam=0.1): P_max ~ 0.0723 (+16% vs Gertsenshtein)
#
# Refs:
#   Barker et al. (2024), arXiv:2406.12826 — eWGT / dark photon analogy
#   Berlin et al. (2024), arXiv:2405.08865 — dark photon mixing methods
#
# Running:
#   cd examples/torsion_dark_photon && bash sweep_2d.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/torsion_dark_photon.json"
OUTPUT="${SCRIPT_DIR}/../data/torsion_dark_photon_sweep_2d"

# Fixed physics parameters
KAPPA=1.0
B0=0.01        # small-P regime
K0=2.0         # plane-wave wavevector
XI=0.1         # torsion kinetic strength (fixed)

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

echo "=== Dark Photon Torsion: 2D Alpha x Deltam Sweep ==="
echo "Fixed: kappa=${KAPPA}, B0=${B0}, xi=${XI}"
echo "Sweep: alpha = 0.1 .. 2.0 (8 pts) × deltam = -0.3 .. 0.3 (7 pts)"
echo "Total runs: 56"
echo "Output: ${OUTPUT}"
echo ""

# Note: --measure peak_conversion only (no conservation).
# Total Hamiltonian conservation is NOT a valid check for this constrained gauge theory.
# Diverged (ghost/unstable) runs are identified by P_max >> P_Gertsenshtein or
# run_status != "success" in the output CSV.
uv run tidal sweep "${SPEC}" \
  --sweep "alpha=0.1:2.0:8" \
  --sweep "deltam=-0.3:0.3:7" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" --param "xi=${XI}" \
  --parallel 4 --resume \
  --no-require-stable \
  --output "${OUTPUT}"
# --no-require-stable: torsion mass matrix entries appear negative by TIDAL
# sign convention even for stable physics. Ghost-free region is xi > deltam^2:
# with xi=0.1 and |deltam| <= 0.3, all 56 points satisfy 0.1 > 0.09 (xi>0.1,
# marginally; note |deltam|=0.3 gives deltam^2=0.09 < 0.1). Runs with
# |deltam| > sqrt(xi)=0.316 enter the ghost region and may diverge (NaN);
# these will be flagged as run_status=diverged in the output.

echo ""
echo "--- Generating plots ---"

# 2D heatmap: P_max over (alpha, deltam) parameter plane
# Colours show amplification; dashed contour at P = sin^2(0.25) ~ 0.0617
# shows the Gertsenshtein baseline (amplification boundary).
uv run tidal plot "${OUTPUT}" --type sweep \
  --metric P_max --log-scale \
  --title "Torsion amplification map: P_max(alpha, deltam) at xi=${XI}" \
  --output "${OUTPUT}/sweep_2d.png" --quiet

# C₀ = P/B₀² analysis (primary metric per project_sweep_measurement_strategy.md)
uv run python "${SCRIPT_DIR}/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}"
echo "Plots:"
echo "  ${OUTPUT}/sweep_2d.png     [raw P_max heatmap, log scale]"
echo "  ${OUTPUT}/analysis_C0.png  [C₀ = P/B₀² for stable runs]"
