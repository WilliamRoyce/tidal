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


def _read_priors(chain_dir: Path, param_names: list[str]) -> list[dict[str, Any]]:
    """Load priors from inference.json if available, else infer arctan_uniform."""
    inf_json = chain_dir / "inference.json"
    if inf_json.is_file():
        try:
            with inf_json.open() as f:
                meta = json.load(f)
            priors = meta.get("priors")
            if priors:
                return priors
        except (json.JSONDecodeError, OSError) as exc:
            log.warning("could not read priors from %s: %s", inf_json, exc)
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
    ]


def _build_stub(chain_dir: Path, param_names: list[str]) -> StubResult:
    """Build a StubResult pointing at the on-disk chain root."""
    chain_root = chain_dir / "_chains" / "tidal"
    priors = _read_priors(chain_dir, param_names)
    return StubResult(
        method="nested",
        param_names=list(param_names),
        priors=priors,
        metadata={"chain_root": str(chain_root), "priors": priors},
    )


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
    for direction, chain_dir in (("amp", amp_path), ("sup", sup_path)):
        stub = _build_stub(chain_dir, params)
        try:
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
