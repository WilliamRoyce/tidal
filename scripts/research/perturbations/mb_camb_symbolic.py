"""Print CAMB's own scalar perturbation equations in CAMB's variables (R-C, #567).

Why: ``docs/cosmology/conventions.md`` says CAMB's perturbation-variable definitions govern
the seam, and the R-C reproduction of the Ma & Bertschinger equations must be transcribed
into our convention and checked by machine, not by a written table. ``camb.symbolic`` carries
the ΛCDM scalar constraint and evolution equations in covariant notation and can evaluate them
in the Newtonian gauge, the synchronous gauge and the CDM frame. This script prints them
verbatim (``sympy`` ``srepr`` plus a readable form) so ``wolfram/repro_b_mb_scalars.wls`` and ``wolfram/probe_c_signature.wls`` can cite the
printed lines and the kernel can compare xPand's output against them.

Output is deterministic for a given CAMB version; the committed copy
``camb_symbolic_newtonian_2.0.4.txt`` was produced by this script and is regenerated to verify.

Usage: ``uv run python scripts/research/perturbations/mb_camb_symbolic.py [OUT_FILE]``
"""
# cspell:words srepr adotoa dgrho dgpi hdot etak Newt

from __future__ import annotations

import sys
from pathlib import Path

import camb
import sympy
from camb import symbolic as cs


def _section(title: str) -> list[str]:
    return ["", f"## {title}", ""]


def _eqs(label: str, eqs, gauge) -> list[str]:
    lines = [f"### {label}", ""]
    for i, eq in enumerate(eqs):
        res = gauge(eq) if gauge is not None else eq
        if res is True:
            lines.append(f"[{i}] (identically satisfied in this gauge)")
            continue
        if isinstance(res, list):
            res = res[0] if len(res) == 1 else res
        shown = sympy.simplify(res) if isinstance(res, sympy.Eq) else res
        lines.extend([f"[{i}] {shown}", f"    srepr: {sympy.srepr(res)}"])
    lines.append("")
    return lines


def main(out_file: str | None) -> None:
    lines = [
        "# camb.symbolic scalar perturbation equations, printed by mb_camb_symbolic.py",
        "# cspell:words dgrho hdot srepr adotoa etak Newt",
        f"camb {camb.__version__}; sympy {sympy.__version__}",
        "Line element declared by camb.symbolic (newtonian_gauge docstring):",
        "  ds^2 = a^2 ( (1 + 2 Psi_N) dt^2 - (1 - 2 Phi_N) delta_ij dx^i dx^j ),  t = conformal time",
        "Conventions: kappa = 8 pi G; H = adotoa (conformal Hubble); K_fac = 1 - 3K/k^2;",
        "  delta = dgrho/(kappa a^2) = total density perturbation, q = heat flux, Pi = anisotropic stress.",
    ]
    lines += _section("Background")
    lines += _eqs(
        "Friedmann and background substitutions",
        [cs.Friedmann, *cs.Friedmann_subs, *cs.background_eqs],
        None,
    )
    lines += _section("Covariant (frame-independent) forms")
    lines += _eqs("constraints (cons1..cons4 = 0)", cs.constraints, None)
    lines += _eqs("evolution equations (pert_eqs)", cs.pert_eqs, None)
    lines += _eqs("total matter equations (total_eqs)", cs.total_eqs, None)
    lines += _section("Newtonian gauge (sigma = 0; variables Phi_N, Psi_N)")
    lines += _eqs(
        "Newtonian variable definitions (Newt_vars, Newtonian_var_subs)",
        [*cs.Newt_vars, *cs.Newtonian_var_subs],
        None,
    )
    lines += _eqs("constraints", cs.constraints, cs.newtonian_gauge)
    lines += _eqs("evolution equations", cs.pert_eqs, cs.newtonian_gauge)
    lines += _eqs("total matter equations", cs.total_eqs, cs.newtonian_gauge)
    lines += _section("Synchronous gauge (variables eta_s, hdot_s)")
    lines += _eqs("constraints", cs.constraints, cs.synchronous_gauge)
    lines += _eqs("evolution equations", cs.pert_eqs, cs.synchronous_gauge)
    lines += _eqs("total matter equations", cs.total_eqs, cs.synchronous_gauge)
    lines += _section(
        "CDM frame (v_c = 0, A = 0; synchronous gauge in covariant names)"
    )
    lines += _eqs("constraints", cs.constraints, cs.cdm_gauge)
    lines += _eqs("evolution equations", cs.pert_eqs, cs.cdm_gauge)
    text = "\n".join(lines) + "\n"
    if out_file:
        Path(out_file).write_text(text, encoding="utf-8")
        print(f"wrote {out_file} ({len(lines)} lines)")
    else:
        print(text)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
