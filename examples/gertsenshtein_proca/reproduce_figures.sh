#!/usr/bin/env bash
# Reproduces the plasma-Gertsenshtein deliverable figures.
#
# Outputs three publication-quality PNGs demonstrating the v0.30.1
# validation of graviton-photon conversion with an effective photon mass
# against the Raffelt & Stodolsky (1988) two-state mixing formula.
#
#   Figure 1 — 1D Raffelt-Stodolsky Lorentzian at fixed B₀=0.10
#              (from sweep_resonance_1d.sh)
#   Figure 2 — 2D resonance map P_max(B₀, mA²) with predicted
#              resonance line (from sweep_resonance.sh)
#   Figure 3 — mA² family overlay of conversion probability P(t)
#              across on- and off-resonance mA² values
#              (from sweep_mA2_family.sh)
#
# HPC workflow: by default this script REUSES existing sweep output
# directories under examples/data/, so sweeps can run on a cluster and
# the results rsync'd back for plotting. Use --fresh to force re-run.
#
# References
# ----------
#   Gertsenshtein (1962)        JETP 14:84          — original conversion
#   Boccaletti et al. (1970)    Nuovo Cim B 70:129  — localized-field form
#   Raffelt & Stodolsky (1988)  Phys Rev D 37:1237  — two-state mixing
#   Ejlli (2020)                arXiv:2004.02714    — exact B₀ solution
#   Hwang & Noh (2023)          arXiv:2310.04150    — curved-spacetime EM
#   Dandoy et al. (2024)        arXiv:2406.17853    — canonical normalization
#   Domcke et al. (2025)        arXiv:2507.16609    — 3D scattering
#
# Physics parameters
# ------------------
#   kappa = 1, B₀ = 0.10, k = 2π·32/100 = 2.0106, D = 50
#   Coupling:   2·κ·B₀·k = 0.4021
#   Resonance:  mA² = κ²·B₀²/2 = 0.005
#   HWHM:       2·κ·B₀·k = 0.4021  (in Δm² = mA² - 0.005)
#
# Usage
# -----
#   cd examples/gertsenshtein_proca
#   bash reproduce_figures.sh           # reuse existing sweep data
#   bash reproduce_figures.sh --fresh   # force rerun all sweeps
#
# Runtime: ~2 min if data is cached, ~25 min for fresh run

set -euo pipefail
cd "$(dirname "$0")"

FRESH=0
if [[ "${1:-}" == "--fresh" ]]; then
  FRESH=1
fi

DATA_DIR_1D=../data/sweep_proca_resonance_1d
DATA_DIR_2D=../data/sweep_proca_resonance
DATA_DIR_FAMILY=../data/sweep_proca_mA2_family
FIG_DIR=../data/gertsenshtein_proca_figures

mkdir -p "$FIG_DIR"

echo "=== Plasma-Gertsenshtein Figure Reproduction ==="
echo ""
echo "References:"
echo "  Raffelt & Stodolsky 1988, Phys Rev D 37:1237 (two-state mixing)"
echo "  Hwang & Noh 2023, arXiv:2310.04150 (curved-spacetime EM conventions)"
echo "  Domcke et al 2025, arXiv:2507.16609 (3D scattering extension)"
echo ""
echo "Expected sweep data directories:"
echo "  1D scan:       $DATA_DIR_1D"
echo "  2D resonance:  $DATA_DIR_2D"
echo "  mA² family:    $DATA_DIR_FAMILY"
echo ""

# ---------------------------------------------------------------------
# Figure 1 — 1D Raffelt-Stodolsky Lorentzian
# ---------------------------------------------------------------------
echo "=== Figure 1: 1D Raffelt-Stodolsky Lorentzian ==="

if [[ $FRESH -eq 1 || ! -f "$DATA_DIR_1D/results.csv" ]]; then
  if [[ $FRESH -eq 1 ]]; then
    rm -rf "$DATA_DIR_1D"
  fi
  echo "  Running sweep_resonance_1d.sh (40 sims)..."
  bash sweep_resonance_1d.sh
else
  echo "  Reusing existing $DATA_DIR_1D/ (--fresh to force rerun)"
fi

# Lorentzian overlay: P = coupling^2 / (coupling^2 + (mA2 - resonance)^2)
# kappa=1, B0=0.10, k=2.0106 => coupling = 0.4021, resonance = 0.005
tidal plot "$DATA_DIR_1D" --type sweep \
  --metric P_max \
  --title "Plasma Gertsenshtein: P_max(mA^2) at B0=0.10 — Raffelt & Stodolsky 1988" \
  --overlay "(2*kappa*B0*2.0106)**2 / ((2*kappa*B0*2.0106)**2 + (mA2 - kappa**2*B0**2/2)**2)" \
  --output "$FIG_DIR/fig1_lorentzian_1d.png" --quiet

echo "  → $FIG_DIR/fig1_lorentzian_1d.png"
echo "    Physics: peak at mA²=0.005, HWHM at Δm²=2·κ·B₀·k=0.402"
echo ""

# ---------------------------------------------------------------------
# Figure 2 — 2D Resonance Map
# ---------------------------------------------------------------------
echo "=== Figure 2: 2D Resonance Map ==="

if [[ $FRESH -eq 1 || ! -f "$DATA_DIR_2D/results.csv" ]]; then
  if [[ $FRESH -eq 1 ]]; then
    rm -rf "$DATA_DIR_2D"
  fi
  echo "  Running sweep_resonance.sh (300 sims)..."
  bash sweep_resonance.sh
else
  echo "  Reusing existing $DATA_DIR_2D/ (--fresh to force rerun)"
fi

tidal plot "$DATA_DIR_2D" --type sweep \
  --metric P_max \
  --title "Plasma Gertsenshtein: P_max(B0, mA^2) — resonance line mA^2 = kappa^2 B0^2 / 2" \
  --overlay "(2*kappa*B0*2.0106)**2 / ((2*kappa*B0*2.0106)**2 + (mA2 - kappa**2*B0**2/2)**2)" \
  --output "$FIG_DIR/fig2_resonance_map_2d.png" --quiet

echo "  → $FIG_DIR/fig2_resonance_map_2d.png"
echo "    Physics: predicted resonance line mA²=κ²B₀²/2 traces the peak"
echo ""

# ---------------------------------------------------------------------
# Figure 3 — mA² Family Oscillation Traces
# ---------------------------------------------------------------------
echo "=== Figure 3: mA² Family Oscillations ==="

if [[ $FRESH -eq 1 || ! -f "$DATA_DIR_FAMILY/results.csv" ]]; then
  if [[ $FRESH -eq 1 ]]; then
    rm -rf "$DATA_DIR_FAMILY"
  fi
  echo "  Running sweep_mA2_family.sh (6 sims)..."
  bash sweep_mA2_family.sh
else
  echo "  Reusing existing $DATA_DIR_FAMILY/ (--fresh to force rerun)"
fi

# Use sweep-compare to overlay P_max vs mA² with the Lorentzian prediction
tidal plot "$DATA_DIR_FAMILY" --type sweep \
  --metric P_max \
  --title "Plasma Gertsenshtein: mA^2 family at B0=0.10 (transition on→off resonance)" \
  --overlay "(2*kappa*B0*2.0106)**2 / ((2*kappa*B0*2.0106)**2 + (mA2 - kappa**2*B0**2/2)**2)" \
  --output "$FIG_DIR/fig3_mA2_family.png" --quiet

echo "  → $FIG_DIR/fig3_mA2_family.png"
echo "    Physics: 6 discrete mA² points spanning on-resonance to P_max ~ 0.15"
echo ""

# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------
echo "=== Done ==="
echo ""
echo "Deliverable figures in $FIG_DIR/ :"
ls -la "$FIG_DIR"/*.png 2>/dev/null || echo "  (no PNGs produced — check errors above)"
echo ""
echo "Docs: see docs/tex/gertsenshtein_proca.tex for the full derivation and"
echo "      physics motivation, and docs/references.md for the pre-arXiv"
echo "      primary citations (Gertsenshtein 1962, Boccaletti 1970, Raffelt-"
echo "      Stodolsky 1988)."
