#!/usr/bin/env bash
# scripts/listings/render_app_e_symbolic.sh
#
# Driver: regenerate the symbolic-output listings for Appendix E
# (manuscript/sections/appendices/symbolic_outputs.tex).
#
# For each campaign theory listed in the MANIFEST below, this script
# invokes `tidal inspect --latex --latex-format align` on the derived
# JSON spec, strips the wrapping `\begin{align} … \end{align}` lines,
# and writes the bare aligned equation body to
#   manuscript/sections/appendices/listings/eom_<tag>_full.tex
# so that the appendix can `\input` the file inside its own
# `\begin{aligned} … \end{aligned}` figure float.
#
# Idempotent: byte-identical input JSONs produce byte-identical
# output. Planned theories without a derived JSON yet are written
# as one-line placeholders so the appendix's section numbering is
# stable as the campaign progresses.
#
# Usage: bash scripts/listings/render_app_e_symbolic.sh
set -euo pipefail

cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"

OUT_DIR="manuscript/sections/appendices/listings"
mkdir -p "$OUT_DIR"

# Manifest: tag | JSON path | status
# - tag      : output filename stem (eom_<tag>_full.tex)
# - JSON     : derived spec path (relative to repo root); empty if no JSON yet
# - status   : "surveyed" (must have JSON) or "planned" (placeholder ok)
MANIFEST=(
  "baseline|examples/data/gertsenshtein.json|surveyed"
  # Euler--Heisenberg QED 1-loop quartic correction. Non-torsion theory
  # that exercises the perturbative-reduction (Pass 0 / Pass 1) path
  # and reproduces the Adler 7/4 birefringence ratio (cf. validation
  # appendix).
  "euler_heisenberg|examples/data/euler_heisenberg.json|surveyed"
  "reference|examples/data/torsion_gertsenshtein_b5_zero.json|surveyed"
  # `mechanism` reads the EXACT (non-LPS) R²-PGT JSON. The sibling
  # `torsion_gertsenshtein.json` is the LPS-reduced canonical form used
  # by `tidal simulate` and the perturbative-baseline tests; the
  # appendix wants the higher-derivative equations explicit, so it
  # reads `_exact.json` instead. Both JSONs are derived from
  # `examples/torsion_gertsenshtein/theory.toml` and `theory_exact.toml`
  # respectively.
  "mechanism|examples/data/torsion_gertsenshtein_exact.json|surveyed"
  "nonminimal|examples/data/torsion_gertsenshtein_nonminimal.json|surveyed"
  "general_nonminimal|examples/data/torsion_gertsenshtein_general_nonminimal.json|surveyed"
  "parity_odd|examples/data/torsion_gertsenshtein_parity_odd.json|surveyed"
  "dark_photon_plasma|examples/data/dark_photon_plasma.json|surveyed"
)

render_one() {
  local tag="$1"
  local json="$2"
  local status="$3"
  local out="$OUT_DIR/eom_${tag}_full.tex"

  if [[ "$status" == "surveyed" ]]; then
    if [[ ! -f "$json" ]]; then
      echo "ERROR: surveyed theory '$tag' is missing its JSON spec: $json" >&2
      return 1
    fi
    # `gather` format wraps each top-level equation (Lagrangian, EOM,
    # Hamiltonian) in its own `\begin{aligned}…\end{aligned}` block,
    # separates them with `\\[1ex]`, and wraps the whole sequence in
    # `\begin{gather*}…\end{gather*}`. We strip the outer `gather*`
    # wrappers so the file can be \input{} inside a `gathered`
    # environment that the appendix shell controls. The per-equation
    # `aligned` structure breaks the global `&=` column drag — each
    # equation centres on its own width, eliminating the "ragged
    # right" whitespace band when short EOMs share a listing with a
    # wide Hamiltonian. The depth-aware wrap stage at --width 400
    # still breaks long lines at top-level term boundaries; the
    # appendix figure float additionally wraps each body in
    # \adjustbox{max width=\linewidth} as a safety net for any
    # residual overflows.
    uv run tidal inspect "$json" --latex --latex-format gather \
        --symbols manuscript/latex_symbols.toml \
      | sed -e '/^%/d' -e '/^\\begin{gather\*}/d' -e '/^\\end{gather\*}/d' \
      | python3 scripts/listings/wrap_long_lines.py --width 400 \
      > "$out"

    # Second pass: emit the kinetic-matrix listing (issue #372). The
    # matrix is assembled from the same JSON via
    # tidal/symbolic/kinetic_matrix.py and rendered as a bmatrix by
    # tidal/symbolic/latex.py::kinetic_matrix_to_latex. Each theory
    # gets eom_<tag>_kinetic_matrix.tex alongside the EOM listing.
    local km_out="$OUT_DIR/eom_${tag}_kinetic_matrix.tex"
    uv run tidal inspect "$json" --latex --latex-format kinetic-matrix \
        --symbols manuscript/latex_symbols.toml \
      > "$km_out"
  else
    # Placeholder for a theory whose derivation has not yet been run.
    # Keeps the subsection's \input target valid; replaced when the
    # JSON lands.
    printf '%%%% Placeholder: %s — derivation pending.\n\\text{Equations of motion to be added when the derivation completes.}\n' "$tag" > "$out"
    printf '%%%% Placeholder: %s kinetic matrix — derivation pending.\n\\text{Kinetic matrix to be added when the derivation completes.}\n' "$tag" > "$OUT_DIR/eom_${tag}_kinetic_matrix.tex"
  fi
  local lines
  lines=$(wc -l <"$out")
  printf '  wrote %-40s (%s lines)\n' "$out" "$lines"
}

echo "Rendering symbolic-output listings for Appendix E"
for entry in "${MANIFEST[@]}"; do
  IFS='|' read -r tag json status <<<"$entry"
  render_one "$tag" "$json" "$status"
done

# Post-render notation patch: the renderer emits the graviton metric
# perturbation as \mathcal{H} and the photon vector potential as
# \mathcal{A} (the default `_calligraphic_head` behaviour in
# tidal/symbolic/latex.py). The appendix convention is lowercase italic
# h_{\mu\nu} (standard GR / TT-gauge notation) and italic capital A_\mu
# (standard EM notation), with torsion left calligraphic to keep it
# visually distinct. We re-letter the listings in place rather than
# changing the renderer so that other consumers (tests, ad-hoc inspect
# calls) keep the historical default.
sed -i \
  -e 's/\\mathcal{H}/h/g' \
  -e 's/\\mathcal{A}/A/g' \
  "$OUT_DIR"/eom_*_full.tex \
  "$OUT_DIR"/eom_*_kinetic_matrix.tex

# Refresh the caption-count sidecar (EOM / constraint / term counts per
# theory). Idempotent; reads only the listings written above.
python3 scripts/listings/compute_appendix_counts.py

echo "Done. Listings live in $OUT_DIR/."
