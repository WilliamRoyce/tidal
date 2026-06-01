r"""Emit four TeX fragments summarizing the §4 D_KL parameter-importance
analysis, sourced from the parameter_importance.json files written by
`recompute_parameter_kl.py`.

Outputs (written to `manuscript/sections/listings/`):
  d_kl_summary.tex          — survey-level table (one row per theory)
  d_kl_per_theory.tex       — per-theory ranked sub-tables
  d_kl_headline_claims.tex  — one-sentence headline claim per theory
  d_kl_methods_appendix.tex — methodology paragraph for App J

Each fragment is self-contained TeX that `\input`s cleanly into the manuscript.
The headline-claims fragment is the starter set for the §4 per-theory
sentences; final manuscript wording can edit the prose around the numeric
values.

Run from the repo root:
    uv run python scripts/analysis/extract_d_kl.py
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

try:
    import tomllib  # type: ignore[import]
except ImportError:
    import tomli as tomllib  # type: ignore[import]

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("extract_d_kl")


# ----- Parameter name → display label (LaTeX math mode) ----------------------

PARAM_LABEL: dict[str, str] = {
    # Torsion-squared
    "beta1": r"\beta_1",
    "beta2": r"\beta_2",
    "beta3": r"\beta_3",
    "alpha1": r"\alpha_1",
    "alpha2": r"\alpha_2",
    "alpha3": r"\alpha_3",
    # Barker single χ + ξ
    "xi": r"\xi",
    "chi": r"\chi",
    "delta1": r"\delta_1",
    # Shapiro ζ
    "zeta1": r"\zeta_1",
    "zeta2": r"\zeta_2",
    "zeta3": r"\zeta_3",
    # Parity-odd minimal
    "d14": r"d_{14}",
    "d15": r"d_{15}",
    "d17": r"d_{17}",
    "d19": r"d_{19}",
    "d20": r"d_{20}",
    "d21": r"d_{21}",
    # Dark-photon plasma legacy names
    "mA2": r"m_A^{2}",
    "deltam": r"\delta_m",
    # EH+Gert
    "rho": r"\rho",
    "sigma": r"\sigma",
}


def _label(name: str) -> str:
    """Return $-wrapped TeX label for a parameter name."""
    if name in PARAM_LABEL:
        return f"${PARAM_LABEL[name]}$"
    # Generated families: chi_i, xi_i, tildechi_i, tildezeta_i
    for prefix, sym in (
        ("tildechi", r"\tilde\chi"),
        ("tildezeta", r"\tilde\zeta"),
        ("chi", r"\chi"),
        ("xi", r"\xi"),
        ("zeta", r"\zeta"),
        ("beta", r"\beta"),
        ("alpha", r"\alpha"),
        ("d", r"d"),
    ):
        if name.startswith(prefix) and name[len(prefix) :].isdigit():
            idx = int(name[len(prefix) :])
            return f"${sym}_{{{idx}}}$"
    return rf"$\mathtt{{{name}}}$"


# ----- Theory name → manuscript display name --------------------------------

THEORY_DISPLAY: dict[str, str] = {
    "dark_photon_plasma": "Dark-photon plasma",
    "ricci_em": "Ricci--EM cross-term",
    "ym_pgt_bahamonde": "Yang--Mills--PGT (Bahamonde)",
    "ym_pgt_barker": "Yang--Mills--PGT (Barker)",
    "ym_pgt_shapiro": "Yang--Mills--PGT (Shapiro)",
    "ym_pgt_full_nonminimal": "Yang--Mills--PGT (full non-minimal)",
    "np_t5": "Non-propagating torsion ($\\xi=0$ control)",
    "chi_closure": "Complete parity-even ($\\chi$-closure)",
    "np_chi_closure": "Non-propagating torsion ($\\xi=0$, $\\chi$-closure)",
    "xi_kinetic_closure": "$\\xi$-kinetic closure",
    "parity_odd_minimal": "Parity-odd minimal",
    "parity_odd_full": "Parity-odd full",
    "complete_parity_odd": "Complete parity-odd",
    "eh_gertsenshtein": "Euler--Heisenberg + Gertsenshtein",
}


def _load_results(config_path: Path) -> list[dict[str, Any]]:
    """Load each theory's parameter_importance.json keyed by theory name."""
    with config_path.open("rb") as f:
        cfg = tomllib.load(f)
    out: list[dict[str, Any]] = []
    for entry in cfg.get("chain", []):
        amp_path = REPO_ROOT / entry["amp_path"]
        json_path = amp_path / "parameter_importance.json"
        if not json_path.is_file():
            log.warning(
                "[%s] no parameter_importance.json yet (run recompute first); skipping",
                entry["theory"],
            )
            continue
        with json_path.open() as f:
            data = json.load(f)
        # Merge config metadata in (geometry, jobid) so the emitter has it.
        data["amp_path"] = entry["amp_path"]
        data["sup_path"] = entry["sup_path"]
        out.append(data)
    return out


def _fmt(x: float, decimals: int = 2) -> str:
    """Format a float with NaN tolerance."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    return f"{x:.{decimals}f}"


def _top_k(d: dict[str, float], k: int) -> list[tuple[str, float]]:
    return sorted(
        ((n, v) for n, v in d.items() if math.isfinite(v)),
        key=lambda kv: -kv[1],
    )[:k]


def emit_summary(results: list[dict[str, Any]], out_path: Path) -> None:
    """Emit the survey-level D_KL summary table.

    Columns: theory, dim, logZ_+, logZ_-, D_KL_+, D_KL_-, d_G_+.
    """
    lines: list[str] = []
    lines.extend(
        (
            r"% Auto-generated by scripts/analysis/extract_d_kl.py — do not hand-edit.",
            r"\begin{tabular}{lccccccc}",
            r"\hline\hline",
            r"Theory & $d$ & $\log\mathcal{Z}_+ \pm \sigma$ & $\log\mathcal{Z}_- \pm \sigma$ & $D_\mathrm{KL,+}$ & $D_\mathrm{KL,-}$ & $d_{G,+}$ & top-$D_\mathrm{KL,\times}$ \\",
            r"\hline",
        )
    )
    for r in results:
        theory = THEORY_DISPLAY.get(r["theory"], r["theory"])
        amp = r.get("amp", {})
        sup = r.get("sup", {})
        cross_top = _top_k(r.get("cross_amp_sup_kl", {}), 1)
        cross_str = (
            f"{_label(cross_top[0][0])} ({_fmt(cross_top[0][1])})"
            if cross_top
            else "--"
        )
        lines.append(
            f"{theory} & {r['n_params']}"
            f" & ${_fmt(amp.get('log_evidence'))} \\pm {_fmt(amp.get('log_evidence_err'))}$"
            f" & ${_fmt(sup.get('log_evidence'))} \\pm {_fmt(sup.get('log_evidence_err'))}$"
            f" & ${_fmt(amp.get('d_kl'))} \\pm {_fmt(amp.get('d_kl_err'))}$"
            f" & ${_fmt(sup.get('d_kl'))} \\pm {_fmt(sup.get('d_kl_err'))}$"
            f" & ${_fmt(amp.get('d_g'))}$"
            f" & {cross_str} \\\\"
        )
    lines.extend((r"\hline\hline", r"\end{tabular}"))
    out_path.write_text("\n".join(lines) + "\n")
    log.info("wrote %s (%d theories)", out_path.relative_to(REPO_ROOT), len(results))


def emit_per_theory(results: list[dict[str, Any]], out_path: Path) -> None:
    r"""Per-theory ranked sub-tables.

    Top-$k$ parameters by amp marginal $D_\\mathrm{KL}$ with side-by-side
    amp / sup / cross columns.
    """
    lines: list[str] = []
    lines.append(
        r"% Auto-generated by scripts/analysis/extract_d_kl.py — do not hand-edit."
    )
    for r in results:
        theory = THEORY_DISPLAY.get(r["theory"], r["theory"])
        amp_marg = r.get("amp", {}).get("marginal_d_kl", {})
        sup_marg = r.get("sup", {}).get("marginal_d_kl", {})
        cross = r.get("cross_amp_sup_kl", {})
        # Rank by union score: max of (amp marginal, sup marginal, cross)
        names = list(amp_marg.keys())
        scores = {
            n: max(
                amp_marg.get(n, 0.0) or 0.0,
                sup_marg.get(n, 0.0) or 0.0,
                cross.get(n, 0.0) or 0.0,
            )
            for n in names
        }
        k = 10 if r["n_params"] > 10 else min(5, r["n_params"])
        ranked = sorted(names, key=lambda n: -scores.get(n, 0.0))[:k]

        lines.extend(
            (
                rf"\paragraph*{{{theory}}}",
                r"\begin{tabular}{lccc}",
                r"\hline\hline",
                r"Parameter & $D_\mathrm{KL,+}$ (nats) & $D_\mathrm{KL,-}$ (nats) & $D_\mathrm{KL}(P_+\|P_-)$ \\",
                r"\hline",
            )
        )
        for n in ranked:
            lines.append(
                f"{_label(n)} & {_fmt(amp_marg.get(n))}"
                f" & {_fmt(sup_marg.get(n))} & {_fmt(cross.get(n))} \\\\"
            )
        lines.extend((r"\hline\hline", r"\end{tabular}", ""))
    out_path.write_text("\n".join(lines) + "\n")
    log.info("wrote %s (%d theories)", out_path.relative_to(REPO_ROOT), len(results))


def emit_headline_claims(results: list[dict[str, Any]], out_path: Path) -> None:
    """One-sentence-per-theory factual claim ready to paste into §4."""
    lines: list[str] = []
    lines.append(
        r"% Auto-generated by scripts/analysis/extract_d_kl.py — manually polish wording."
    )
    for r in results:
        theory = THEORY_DISPLAY.get(r["theory"], r["theory"])
        cross = r.get("cross_amp_sup_kl", {})
        amp_marg = r.get("amp", {}).get("marginal_d_kl", {})
        top_cross = _top_k(cross, 2)
        top_amp = _top_k(amp_marg, 2)
        if not top_cross or not top_amp:
            lines.append(f"% {theory}: no carrier data")
            continue

        top_carrier, top_carrier_v = top_cross[0]
        amp_top1, amp_top1_v = top_amp[0]
        cross_str = f"{_label(top_carrier)} ($D_\\mathrm{{KL}}(P_+\\|P_-) = {_fmt(top_carrier_v)}$ nats)"
        amp_str = f"{_label(amp_top1)} ($D_\\mathrm{{KL,+}} = {_fmt(amp_top1_v)}$ nats)"
        next_amp_str = ""
        if len(top_amp) > 1:
            n2, v2 = top_amp[1]
            next_amp_str = f" and {_label(n2)} ($D_\\mathrm{{KL,+}} = {_fmt(v2)}$ nats)"
        # Construct claim sentence.
        claim = (
            rf"% {theory}"
            "\n"
            rf"% Amp carriers (marginal): {amp_str}{next_amp_str}; "
            rf"top cross-amp-sup discriminator: {cross_str}."
        )
        lines.extend((claim, ""))
    out_path.write_text("\n".join(lines) + "\n")
    log.info("wrote %s (%d theories)", out_path.relative_to(REPO_ROOT), len(results))


def emit_methods(out_path: Path) -> None:
    """Static methodology paragraph for App J."""
    text = r"""% Auto-generated by scripts/analysis/extract_d_kl.py — methodology fragment.
\paragraph*{Information gain and per-parameter $D_\mathrm{KL}$}%
\label{KLMethods}
For each surveyed sector we report two flavors of Kullback--Leibler
divergence.  The \emph{per-direction marginal} $D_\mathrm{KL}(P_\pm(\theta_i)
\,\|\, \pi(\theta_i))$ measures, parameter-by-parameter, how much the
amplification ($+$) or suppression ($-$) posterior has contracted relative
to its prior.  It is computed via the standard anesthetic histogram-against-prior
estimator (binning each marginal in the natural space of its prior, where the
prior density is uniform).  Bootstrap uncertainties are propagated from
PolyChord's nested-sampling estimator.

The \emph{cross-amp-sup} divergence $D_\mathrm{KL}(P_+(\theta_i)
\,\|\, P_-(\theta_i))$ measures, parameter-by-parameter, how much the
amplification posterior differs from the suppression posterior.  Both
distributions are estimated as weighted 1D histograms on a shared binning
grid spanning the union of their supports.  This quantity identifies the
parameters that distinguish the two likelihood directions, independent of
the prior baseline.  Together with the directional $\mathrm{ATLAS}$
decomposition of \cref{CubedSphereAtlas}, the marginal and cross divergences
isolate the carrier operators in each surveyed sector; see
\cref{KLPerTheory} for the per-theory tables.
"""
    out_path.write_text(text)
    log.info("wrote %s", out_path.relative_to(REPO_ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO_ROOT / "scripts" / "analysis" / "d_kl_chains.toml",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "manuscript" / "sections" / "listings",
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results = _load_results(args.config)
    if not results:
        log.error("no parameter_importance.json files found")
        return

    emit_summary(results, args.out_dir / "d_kl_summary.tex")
    emit_per_theory(results, args.out_dir / "d_kl_per_theory.tex")
    emit_headline_claims(results, args.out_dir / "d_kl_headline_claims.tex")
    emit_methods(args.out_dir / "d_kl_methods_appendix.tex")


if __name__ == "__main__":
    main()
