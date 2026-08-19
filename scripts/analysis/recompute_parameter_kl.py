r"""Theory-agnostic walker that recomputes parameter_importance for every
chain listed in `d_kl_chains.toml` and writes the result as a sibling
`parameter_importance.json` next to each chain's `inference.json`.

Uses the existing `tidal.inference._importance.compute_parameter_importance`
infrastructure plus the new `compute_cross_kl` to capture:

- Global D_KL ± err (Shannon information gain, posterior vs prior)
- Bayesian model dimensionality d_G ± err
- Per-parameter marginal D_KL (histogram-based against the prior)
- Per-parameter cross-amp-sup D_KL (P_amp vs P_sup)
- logZ ± err

Priors are read from each chain's `inference.json` when available; the
TOML config provides param_names per chain so the column ordering matches
the on-disk chain files.  Add new chains (e.g. localized geometry runs)
by appending [[chain]] entries to `d_kl_chains.toml` — no code edits.

Run from the repo root:
    uv run python scripts/analysis/recompute_parameter_kl.py
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

try:
    import tomllib  # type: ignore[import]
except ImportError:
    import tomli as tomllib  # type: ignore[import]

from tidal.inference._importance import (
    compute_cross_kl,
    compute_parameter_importance,
)

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("recompute_kl")


@dataclass
class StubResult:
    """Minimal duck-typed stand-in for InferenceResult.

    `compute_parameter_importance` only needs `.method`, `.param_names`,
    `.priors` (or `.metadata['priors']`), and `.metadata['chain_root']` so
    `to_anesthetic_samples` can load the chain via `anesthetic.read_chains`.
    """

    method: str
    param_names: list[str]
    metadata: dict[str, Any]
    priors: list[dict[str, Any]]


def _read_priors(
    chain_dir: Path,
    param_names: list[str],
) -> tuple[list[dict[str, Any]], str]:
    """Load priors from inference.json if available, else fabricate arctan.

    Returns ``(records, provenance)`` where ``provenance`` is ``"recorded"``
    or ``"fabricated"``.  Fabrication is LOUD: the Barker amp chain
    (GH #434) had an ``inference.json`` *without* a ``priors`` key and was
    silently scored against fabricated all-arctan(±89) priors for months —
    file existence is not prior provenance, so the caller records the
    provenance in the output JSON and the report prints it from data.
    """
    inf_json = chain_dir / "inference.json"
    if inf_json.is_file():
        try:
            with inf_json.open() as f:
                meta = json.load(f)
            priors = meta.get("priors")
            if priors:
                return priors, "recorded"
            log.warning(
                "[%s] inference.json exists but has NO priors key "
                "(pre-schema or reconstructed run) — fabricating "
                "all-arctan(+-89) priors; add prior_overrides in "
                "d_kl_chains.toml for any parameter whose sbatch template "
                "used a different kind (GH #434)",
                chain_dir.name,
            )
        except (json.JSONDecodeError, OSError) as exc:
            log.warning("could not read priors from %s: %s", inf_json, exc)
    else:
        log.warning(
            "[%s] no inference.json — fabricating all-arctan(+-89) priors; "
            "verify against the sbatch template and add prior_overrides "
            "where wrong (GH #434)",
            chain_dir.name,
        )
    # Fallback: assume arctan_uniform with the convention used elsewhere
    # in the survey (±89 for sign-symmetric couplings).  This is only used
    # to support older chains pre-dating the inference.json schema.
    return [
        {
            "kind": "scalar",
            "name": name,
            "distribution": "arctan_uniform",
            "low": -89.0,
            "high": 89.0,
        }
        for name in param_names
    ], "fabricated"


def _align_priors_to_params(
    priors: list[dict[str, Any]],
    param_names: list[str],
    chain_label: str,
) -> list[dict[str, Any]]:
    """Rename prior records positionally onto the TOML ``param_names``.

    ``d_kl_chains.toml`` may relabel chain columns (e.g. ricci_em's
    ``alpha*`` → ``beta*`` so figure axes match §4 notation) while the
    recorded prior names keep the original labels.  The name-based lookup
    in ``compute_parameter_importance`` then misses every renamed column
    and silently falls back to a uniform-over-sample-range reference
    (GH #420 — the miss was masked pre-fix because the saturated
    estimator produced near-identical numbers either way).

    Records are matched to columns positionally: record order equals
    column order by construction (``build_prior_transform`` concatenates
    prior dimensions in declaration order, and the chain writer emits
    columns the same way).  ``radial_angular`` records consume
    ``len(names)`` consecutive columns.  Any arity mismatch RAISES —
    proceeding with misaligned (or original-named) records would score
    every column against the wrong or no prior, the silent-misassignment
    class the #420 campaign was about; the caller skips the chain with a
    visible error instead.
    """
    aligned: list[dict[str, Any]] = []
    pos = 0
    for rec in priors:
        if rec.get("kind") == "radial_angular":
            names = list(rec.get("names", []))
            new_names = param_names[pos : pos + len(names)]
            if len(new_names) != len(names):
                msg = (
                    f"[{chain_label}] prior/param arity mismatch: joint "
                    f"record needs {len(names)} columns, {len(new_names)} "
                    f"remain — refusing to score against misaligned priors"
                )
                raise ValueError(msg)
            if names != new_names:
                log.info(
                    "[%s] positional prior match: %s -> %s",
                    chain_label,
                    names,
                    new_names,
                )
            aligned.append({**rec, "names": new_names})
            pos += len(names)
        else:
            if pos >= len(param_names):
                msg = (
                    f"[{chain_label}] more prior records than columns — "
                    f"refusing to score against misaligned priors"
                )
                raise ValueError(msg)
            new_name = param_names[pos]
            if rec.get("name") != new_name:
                log.info(
                    "[%s] positional prior match: %r -> %r (%s)",
                    chain_label,
                    rec.get("name"),
                    new_name,
                    rec.get("distribution"),
                )
            aligned.append({**rec, "name": new_name})
            pos += 1
    if pos != len(param_names):
        msg = (
            f"[{chain_label}] {pos} prior dimensions for "
            f"{len(param_names)} columns — refusing to score against "
            f"misaligned priors"
        )
        raise ValueError(msg)
    return aligned


def _apply_prior_overrides(
    priors: list[dict[str, Any]],
    overrides: dict[str, str],
    chain_label: str,
) -> list[dict[str, Any]]:
    """Replace named scalar prior records with ``DIST:LO:HI`` overrides.

    ``d_kl_chains.toml`` entries may carry ``prior_overrides = {name =
    "log_uniform:1e-3:1e3", ...}`` for chains that pre-date the
    inference.json priors schema, where ``_read_priors`` fabricates an
    all-arctan fallback.  The override strings use the same syntax as the
    campaign ``--prior`` flags (see scripts/hpc_templates/*.sbatch, the
    provenance for these values), and are validated through
    :func:`tidal.inference._prior.parse_prior` so a malformed value fails
    with the parameter name and the offending string instead of a bare
    unpacking error.  Names refer to the TOML ``param_names`` (i.e.
    post-alignment).

    Raises
    ------
    ValueError
        If an override string is malformed, or names a parameter that has
        no scalar record (a typo, or a coupling inside a radial_angular
        joint record — which cannot be overridden per-component).
    """
    from tidal.inference._prior import parse_prior

    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    for rec in priors:
        name = rec.get("name")
        if name in overrides:
            try:
                parsed = parse_prior(f"{name}={overrides[name]}")
            except ValueError as exc:
                msg = (
                    f"[{chain_label}] malformed prior_overrides entry "
                    f"{name} = {overrides[name]!r}: {exc}"
                )
                raise ValueError(msg) from exc
            log.info(
                "[%s] prior override: %s = %s (was %s)",
                chain_label,
                name,
                overrides[name],
                rec.get("distribution"),
            )
            out.append(
                {
                    "kind": "scalar",
                    "name": name,
                    "distribution": parsed.distribution,
                    "low": parsed.low,
                    "high": parsed.high,
                }
            )
            seen.add(str(name))
        else:
            out.append(rec)
    missing = set(overrides) - seen
    if missing:
        # An unmatched override means the numbers would silently be wrong
        # for exactly the parameter someone bothered to override — refuse.
        msg = (
            f"[{chain_label}] prior_overrides for parameter(s) with no "
            f"scalar record: {sorted(missing)} (typo, or the coupling "
            f"lives inside a radial_angular joint record)"
        )
        raise ValueError(msg)
    return out


def _build_stub(
    chain_dir: Path,
    param_names: list[str],
    prior_overrides: dict[str, str] | None = None,
) -> tuple[StubResult, str]:
    """Build a StubResult pointing at the on-disk chain root.

    Returns ``(stub, provenance)`` — provenance is ``"recorded"`` /
    ``"fabricated"``, with ``"+overrides"`` appended when
    ``prior_overrides`` replaced any record, and is persisted into the
    output JSON so the report can print prior provenance from data
    instead of inferring it from file existence (GH #434).
    """
    chain_root = chain_dir / "_chains" / "tidal"
    raw_priors, provenance = _read_priors(chain_dir, param_names)
    priors = _align_priors_to_params(raw_priors, param_names, chain_dir.name)
    if prior_overrides:
        priors = _apply_prior_overrides(priors, prior_overrides, chain_dir.name)
        provenance += "+overrides"
    return StubResult(
        method="nested",
        param_names=list(param_names),
        priors=priors,
        metadata={"chain_root": str(chain_root), "priors": priors},
    ), provenance


def _process_chain(entry: dict[str, Any], n_bootstrap: int) -> None:
    theory = entry["theory"]
    amp_path = REPO_ROOT / entry["amp_path"]
    sup_path = REPO_ROOT / entry["sup_path"]
    params = list(entry["param_names"])

    if not (amp_path / "_chains").is_dir() or not (sup_path / "_chains").is_dir():
        log.warning("[%s] skipping — chain dirs missing", theory)
        return
    # PolyChord runs that produced zero dead points (e.g. catastrophic
    # eigenvalue overflow) write only `_chains/tidal_phys_live.txt`, never
    # `_chains/tidal.txt`.  These can't be analyzed via anesthetic; skip.
    if not (amp_path / "_chains" / "tidal.txt").is_file():
        log.warning("[%s] skipping — no dead points (chain did not converge)", theory)
        return

    log.info("[%s] processing %d params", theory, len(params))

    output: dict[str, Any] = {
        "theory": theory,
        "geometry": entry.get("geometry", "unknown"),
        "jobid": entry.get("jobid"),
        "n_params": len(params),
        "param_names": params,
    }

    # Per-direction parameter importance (amp + sup).
    prior_overrides: dict[str, str] = dict(entry.get("prior_overrides", {}))
    for direction, chain_dir in (("amp", amp_path), ("sup", sup_path)):
        try:
            stub, provenance = _build_stub(chain_dir, params, prior_overrides)
            imp = compute_parameter_importance(stub, n_bootstrap=n_bootstrap)
        except (ValueError, RuntimeError, KeyError) as exc:
            log.warning(
                "[%s/%s] compute_parameter_importance failed: %s",
                theory,
                direction,
                exc,
            )
            continue
        output[direction] = {
            "d_kl": imp.d_kl,
            "d_kl_err": imp.d_kl_err,
            "d_g": imp.d_g,
            "d_g_err": imp.d_g_err,
            "log_evidence": imp.log_evidence,
            "log_evidence_err": imp.log_evidence_err,
            "marginal_d_kl": imp.marginal_d_kl,
            # Self-check block (#420): superadditivity, saturation, n_eff.
            "consistency": imp.consistency,
            # recorded | fabricated [+overrides] — from _read_priors, so the
            # report prints provenance from data, not file existence (#434).
            "priors_provenance": provenance,
        }

    # Cross-amp-sup per-parameter D_KL: load both chains via anesthetic
    # directly so we can pass them to `compute_cross_kl`.
    try:
        from anesthetic import read_chains

        amp_ns = read_chains(
            str(amp_path / "_chains" / "tidal"),
            columns=params,
            labels=params,
        )
        sup_ns = read_chains(
            str(sup_path / "_chains" / "tidal"),
            columns=params,
            labels=params,
        )
        cross = compute_cross_kl(amp_ns, sup_ns, params)
        output["cross_amp_sup_kl"] = cross
    except (ImportError, FileNotFoundError, OSError, ValueError) as exc:
        log.warning("[%s] cross-KL failed: %s", theory, exc)
        output["cross_amp_sup_kl"] = {name: float("nan") for name in params}

    # Write the parameter_importance JSON into the amp chain dir
    # (the chain anchor everywhere else; sup is mirrored implicitly).
    out_path = amp_path / "parameter_importance.json"
    with out_path.open("w") as f:
        json.dump(output, f, indent=2, allow_nan=True)
    log.info(
        "[%s] wrote %s (amp D_KL=%.2f, sup D_KL=%.2f)",
        theory,
        out_path.relative_to(REPO_ROOT),
        output.get("amp", {}).get("d_kl", float("nan")),
        output.get("sup", {}).get("d_kl", float("nan")),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO_ROOT / "scripts" / "analysis" / "d_kl_chains.toml",
        help="Path to chain inventory TOML",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=100,
        help="anesthetic stats bootstrap samples (default 100)",
    )
    parser.add_argument(
        "--theory",
        type=str,
        default=None,
        help="If set, process only the chain with this theory name",
    )
    args = parser.parse_args()

    with args.config.open("rb") as f:
        cfg = tomllib.load(f)

    for entry in cfg.get("chain", []):
        if args.theory and entry["theory"] != args.theory:
            continue
        try:
            _process_chain(entry, n_bootstrap=args.n_bootstrap)
        except Exception:
            log.exception("[%s] uncaught error", entry.get("theory", "?"))


if __name__ == "__main__":
    main()
