# Research artifact -- H6 spectrum-module design study (2026-08-30 / 2026-09-03).
# Preserved by H8 from a temporary directory; see this directory's README.md.
# Findings from running it are written up in docs/cosmology/spectrum_design.md
# section 6.1. Paths in the usage note below are historical and no longer exist;
# the production reader lives in tidalcosmo/spectrum/ (see
# docs/cosmology/stage1_engineering_plan.md section 5).
"""
Extract per-sector M0/M1/M2 polynomial-matrix structure from a PSALTer WXF.

Step 1 of the Python pipeline. The Wolfram-side `JuliaExport.m` does the
equivalent today via regex string substitution; this module does it as an
AST translation from the WXF tree into SymPy expressions, kernel-free.

Output is a `Theory` dataclass with:
- theory_name
- couplings: sorted list of sympy.Symbol (matching JuliaExport's c[i] order)
- sectors: list of `Sector` (one per non-empty spin sector)
    - spin_label
    - dim
    - M0, M1, M2 as sympy.Matrix in the couplings (Def-coefficient extracted)

Usage:
    ./jax_test/.venv/bin/python jax_test/wxf_extract.py [path-to-wxf]
"""

from dataclasses import dataclass, field
from pathlib import Path

import sympy as sp
from wolframclient.deserializers import binary_deserialize
from wolframclient.language.expression import WLFunction, WLSymbol

DEF_SYMBOL_NAME = "xAct`PSALTer`Def"

WL_TO_SYMPY_HEAD = {
    "Plus": sp.Add,
    "Times": sp.Mul,
    "Power": sp.Pow,
    "Rational": sp.Rational,
    "Complex": lambda re, im: sp.Add(re, sp.Mul(im, sp.I), evaluate=False),
}


@dataclass
class Sector:
    spin_label: str
    dim: int
    M0: sp.Matrix
    M1: sp.Matrix
    M2: sp.Matrix
    # Degree of det M(k) as polynomial in z = k^2. Equal to dim when M2 is
    # full-rank, smaller when M2 is rank-deficient. -1 if det M(k) is
    # identically zero (gauge-symmetric sector, out of scope here).
    # The wave operator is assumed to be Hermitian for real k (i.e. M0, M2
    # real symmetric and M1 imaginary antisymmetric), so det M(k) is even
    # in k and the polynomial in k^2 is well-defined; this is asserted at
    # extract time and an error is raised if odd-power-in-k coefficients
    # appear (would indicate a non-Hermitian wave operator).
    z_degree: int = -1


@dataclass
class Theory:
    theory_name: str
    couplings: list[sp.Symbol]
    sectors: list[Sector] = field(default_factory=list)


def _short_symbol_name(wl_name: str) -> str:
    """Strip Wolfram context prefixes (e.g. 'Global`Theta1' -> 'Theta1').

    Keep `xAct`PSALTer`Def` intact so we can spot it as the momentum variable.
    """
    if wl_name == DEF_SYMBOL_NAME:
        return wl_name
    return wl_name.rsplit("`", 1)[-1]


def wl_to_sympy(node):
    """Translate a wolframclient WXF node into a sympy expression."""
    if isinstance(node, (int, float)):
        return sp.Integer(node) if isinstance(node, int) else sp.Float(node)
    if isinstance(node, WLSymbol):
        return sp.Symbol(_short_symbol_name(node.name))
    if isinstance(node, WLFunction):
        head_name = node.head.name if isinstance(node.head, WLSymbol) else None
        ctor = WL_TO_SYMPY_HEAD.get(head_name)
        if ctor is None:
            msg = f"Don't know how to translate WLFunction[{head_name}] (n_args={len(node.args)})"
            raise NotImplementedError(msg)
        return ctor(*[wl_to_sympy(a) for a in node.args])
    msg = f"Unhandled node type: {type(node).__name__}"
    raise NotImplementedError(msg)


def matrix_from_wl(rows_tuple) -> sp.Matrix:
    """Translate a WL matrix (tuple of tuples) into a sympy Matrix."""
    return sp.Matrix([[wl_to_sympy(c) for c in row] for row in rows_tuple])


def split_def_powers(M: sp.Matrix, def_sym: sp.Symbol):
    """Decompose M(Def) = M0 + Def*M1 + Def^2*M2.

    Asserts no higher powers exist (PSALTer wave operators are quadratic in
    momentum); raises if that assumption fails.
    """
    M_poly = sp.expand(M.applyfunc(sp.expand))
    M0 = M_poly.applyfunc(lambda e: sp.Poly(e, def_sym).nth(0) if e.has(def_sym) else e)
    M1 = M_poly.applyfunc(
        lambda e: sp.Poly(e, def_sym).nth(1) if e.has(def_sym) else sp.S.Zero
    )
    M2 = M_poly.applyfunc(
        lambda e: sp.Poly(e, def_sym).nth(2) if e.has(def_sym) else sp.S.Zero
    )
    # Assert nothing higher than k^2.
    for i in range(M_poly.rows):
        for j in range(M_poly.cols):
            e = M_poly[i, j]
            if e.has(def_sym):
                degree = sp.Poly(e, def_sym).degree()
                if degree > 2:
                    msg = f"Matrix entry [{i},{j}] has Def^{degree}; expected <= 2"
                    raise AssertionError(msg)
    return M0, M1, M2


SPIN_LABELS = ["0+", "1-", "2+"]


def _z_degree(M0: sp.Matrix, M1: sp.Matrix, M2: sp.Matrix) -> int:
    """Degree of det M(k) as polynomial in z = k^2.

    Asserts the polynomial is purely even in k (i.e. all odd-power
    coefficients vanish symbolically). This holds whenever the wave
    operator is Hermitian for real k, which our pipeline assumes
    throughout (M0, M2 real symmetric and M1 imaginary antisymmetric).

    Returns
    -------
        -1 if det M(k) is identically zero (gauge-symmetric, out of scope).
         0 if det is a non-zero constant in k (no propagating mode).
         d if det is a degree-2d polynomial in k (= degree d in z).
    """
    k = sp.Symbol("__k_z_degree_probe__")
    M_of_k = M0 + k * M1 + (k**2) * M2
    det_expr = M_of_k.det()
    if det_expr == 0:
        return -1
    poly = sp.Poly(det_expr, k)
    if poly.is_zero:
        return -1
    deg_k = poly.degree()
    # Enforce evenness in k.
    for i in range(1, deg_k + 1, 2):
        if poly.nth(i) != 0:
            msg = (
                f"det M(k) has non-zero odd-power coefficient at k^{i}: "
                f"{poly.nth(i)}. The pipeline assumes M(k) is Hermitian "
                f"for real k (M0, M2 real symmetric, M1 imaginary antisymmetric), "
                f"which would make det even in k. Got det = {det_expr}"
            )
            raise ValueError(msg)
    return deg_k // 2


def extract(wxf_path: Path) -> Theory:
    raw = wxf_path.read_bytes()
    obj = binary_deserialize(raw)

    # Top-level Association is keyed by WLSymbol("xAct`PSALTer`<Slot>"), not str.
    slots = {k.name: v for k, v in obj.items()}
    wave_operator_full = slots["xAct`PSALTer`WaveOperator"]

    def_sym = sp.Symbol(DEF_SYMBOL_NAME)

    # Theory name = filename stem with prefix stripped.
    stem = wxf_path.stem
    theory_name = stem.removeprefix("ParticleSpectrograph")

    sectors: list[Sector] = []
    all_couplings: set = set()

    for sector_idx, sector_matrix in enumerate(wave_operator_full):
        if not sector_matrix or not sector_matrix[0]:
            continue
        M = matrix_from_wl(sector_matrix)
        M0, M1, M2 = split_def_powers(M, def_sym)
        for sub in (M0, M1, M2):
            for entry in sub:
                all_couplings.update(entry.free_symbols)
        zd = _z_degree(M0, M1, M2)
        sectors.append(
            Sector(
                spin_label=SPIN_LABELS[sector_idx]
                if sector_idx < len(SPIN_LABELS)
                else f"sector_{sector_idx}",
                dim=M.rows,
                M0=M0,
                M1=M1,
                M2=M2,
                z_degree=zd,
            )
        )

    all_couplings.discard(def_sym)
    couplings = sorted(all_couplings, key=lambda s: s.name)

    return Theory(theory_name=theory_name, couplings=couplings, sectors=sectors)
