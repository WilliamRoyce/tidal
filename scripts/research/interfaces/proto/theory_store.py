"""Shared theory loading, fingerprinting and derived-spectrum store for the R-1 prototypes.

Research prototype for R-1 (#566), not production code. It exercises the design recorded
in ``docs/cosmology/r1_planning_record.md`` §3.1-§3.3:

* the theory definition (``fields``, ``lagrangian``, ``derivation``) is read from YAML with
  **Cobaya's own file loader**, so ``!defaults`` includes expand exactly as ``cobaya-run``
  expands them;
* the fingerprint hashes only the theory-defining keys, canonically serialized, plus a
  schema version — never priors, likelihoods or sampler settings;
* derived spectra live content-addressed in a *data* store (``$XDG_DATA_HOME`` semantics,
  not cache: they cannot be regenerated without Wolfram), searched along an env-var path
  whose first entry is writable, with a manifest and a sha256 verified on read.

Nothing here starts a Wolfram kernel.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

from cobaya.yaml import yaml_load_file

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCHEMA = "r1-proto-1"
"""Fingerprint recipe version. Bump on any change to what is hashed or how."""

THEORY_CLASS = "tidalcosmo.SpectatorTheory"
"""The component name a run file uses; the prototype gate registers under its own name."""

THEORY_KEYS = ("fields", "lagrangian", "derivation")
DERIVATION_DEFAULTS = {"max_laurent_depth": 1}

CONVENTIONS = {"signature": [1, -1, -1, -1], "epsilon0123": 1}
"""PSALTer's conventions (PSALTer.m:101, DefGeometry.m:69-72); asserted on every read."""

STORE_ENV = "TIDALCOSMO_SPECTRA_PATH"


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def theory_from_mapping(block: dict[str, Any]) -> dict[str, Any]:
    """Return the canonical theory definition from a theory mapping.

    Structural checks only (types and keys). Symbolic validation of operator strings is the
    Wolfram package's job and never happens in Python.
    """
    unknown = set(block) - set(THEORY_KEYS)
    if unknown:
        msg = f"unknown theory keys {sorted(unknown)}; allowed: {list(THEORY_KEYS)}"
        raise ValueError(msg)
    fields = block.get("fields")
    lagrangian = block.get("lagrangian")
    if not isinstance(fields, dict) or not fields:
        raise ValueError("theory needs a non-empty 'fields' mapping")
    if not isinstance(lagrangian, dict) or not lagrangian:
        msg = "theory needs a non-empty 'lagrangian' mapping of coupling -> operator string"
        raise ValueError(msg)
    for coupling, operator in lagrangian.items():
        if not isinstance(coupling, str) or not isinstance(operator, str):
            msg = (
                f"lagrangian entry {coupling!r}: key and operator must both be strings"
            )
            raise TypeError(msg)
    derivation = {**DERIVATION_DEFAULTS, **(block.get("derivation") or {})}
    return {"fields": fields, "lagrangian": lagrangian, "derivation": derivation}


def load_theory(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Load a theory definition from a theory file or from a Cobaya run file.

    A run file carries it under ``theory: <component>: model:`` (usually via ``!defaults``).
    """
    info = yaml_load_file(str(path))
    if "theory" in info:
        components = info["theory"] or {}
        blocks = [
            b for b in components.values() if isinstance(b, dict) and "model" in b
        ]
        if len(blocks) != 1:
            msg = f"{path}: expected exactly one theory component with a 'model:' block"
            raise ValueError(msg)
        return theory_from_mapping(blocks[0]["model"])
    return theory_from_mapping(info)


# ---------------------------------------------------------------------------
# Fingerprint
# ---------------------------------------------------------------------------


def canonical_json(theory: dict[str, Any]) -> str:
    """Serialize the theory-defining content deterministically (sorted keys, no spaces)."""
    payload = {"schema": SCHEMA, **theory}
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    )


def fingerprint(theory: dict[str, Any]) -> str:
    """sha256 of the canonical theory. The user never sees or types this."""
    return hashlib.sha256(canonical_json(theory).encode("utf-8")).hexdigest()


def couplings(theory: dict[str, Any]) -> list[str]:
    """The coupling roster, in a stable order."""
    return sorted(theory["lagrangian"])


# ---------------------------------------------------------------------------
# Store
# ---------------------------------------------------------------------------


def default_store() -> Path:
    """User data directory, per XDG: data, not cache (see the planning record §3.3)."""
    base = os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local" / "share")
    return Path(base) / "tidalcosmo" / "spectra"


def search_path(spectra_path: str | None = None) -> list[Path]:
    """Stores in lookup order: explicit option, env-var list, then the default."""
    stores: list[Path] = []
    if spectra_path:
        stores.append(Path(spectra_path))
    stores.extend(Path(p) for p in os.environ.get(STORE_ENV, "").split(os.pathsep) if p)
    stores.append(default_store())
    return stores


def entry_dir(store: Path, fp: str) -> Path:
    return store / f"v{SCHEMA}" / fp[:2] / fp


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_entry(
    fp: str, spectra_path: str | None = None
) -> tuple[Path | None, list[Path]]:
    """First store entry for ``fp`` that has both files, plus the stores searched."""
    stores = search_path(spectra_path)
    for store in stores:
        entry = entry_dir(store, fp)
        if (entry / "manifest.json").is_file() and (entry / "spectrum.wxf").is_file():
            return entry, stores
    return None, stores


def verify_entry(entry: Path, theory: dict[str, Any]) -> dict[str, Any]:
    """Check integrity and identity of a store entry; return its manifest or raise."""
    manifest = json.loads((entry / "manifest.json").read_text())
    problems = []
    if manifest.get("fingerprint") != fingerprint(theory):
        problems.append("fingerprint in manifest does not match this theory")
    if manifest.get("spectrum_sha256") != sha256_file(entry / "spectrum.wxf"):
        problems.append(
            "spectrum.wxf sha256 does not match its manifest (partial or edited file)"
        )
    if manifest.get("conventions") != CONVENTIONS:
        problems.append(
            f"conventions {manifest.get('conventions')} != required {CONVENTIONS}"
        )
    if manifest.get("schema") != SCHEMA:
        problems.append(f"schema {manifest.get('schema')!r} != {SCHEMA!r}")
    if manifest.get("couplings") != couplings(theory):
        problems.append("coupling roster in manifest differs from the theory's")
    if problems:
        raise ValueError("; ".join(problems))
    return manifest
