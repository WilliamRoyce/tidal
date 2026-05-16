#!/usr/bin/env python3
"""Reconstruct inference.json from a truncated PolyChord chain directory.

When a PolyChord INTR job is killed by SLURM before completion, the Python
post-processing code never runs and inference.json is never written.  This
script reads the existing _chains/ files via anesthetic (which handles
truncated dead+live point files correctly) and writes the same outputs that
tidal sample would have produced.

Usage:
    uv run python scripts/reconstruct_polychord_inference.py hpc_results/29232780

Output:
    <output_dir>/inference.json   — summary stats (ESS, MAP, credible intervals)
    <output_dir>/results.csv      — weighted samples
    <output_dir>/results.json     — same data, JSON format

The parameter_importance (D_KL) section is omitted because the prior transform
is not available without the original CLI arguments.  All other fields are
identical to a normally-completed run.  inference.json will include:
    "reconstructed": true
    "n_dead_pts": <count of dead points at truncation>

After running this script, tidal plot --type corner will work on the directory.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def _find_chain_root(output_dir: Path) -> str:
    chains_dir = output_dir / "_chains"
    if not chains_dir.exists():
        sys.exit(f"error: no _chains/ directory found in {output_dir}")
    dead_file = chains_dir / "tidal_dead-birth.txt"
    if not dead_file.exists():
        sys.exit(f"error: {dead_file} not found — nothing to reconstruct")
    return str(chains_dir / "tidal")


def _read_paramnames(chain_root: str) -> list[str] | None:
    """Read parameter names from tidal.paramnames if present."""
    pf = Path(chain_root + ".paramnames")
    if not pf.exists():
        return None
    names = []
    for line in pf.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if parts:
            names.append(parts[0])
    return names or None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reconstruct inference.json from a truncated PolyChord chain."
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing _chains/ (e.g. hpc_results/29232780)",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Optional label appended to inference.json metadata",
    )
    parser.add_argument(
        "--param-names",
        default=None,
        metavar="NAME[,NAME,...]",
        dest="param_names",
        help=(
            "Comma-separated parameter names (e.g. 'beta1,beta2,xi,delta1'). "
            "Required when tidal.paramnames is absent (truncated chain without "
            "write_paramnames output). Count must match ndim in the chain."
        ),
    )
    args = parser.parse_args()

    output_path = Path(args.output_dir).resolve()
    chain_root = _find_chain_root(output_path)

    try:
        from anesthetic import read_chains
    except ImportError:
        sys.exit("error: anesthetic is not installed. Run: uv sync")

    print(f"Reading chains from: {chain_root}_dead-birth.txt")
    ns = read_chains(chain_root)

    # Determine parameter names and ndim.
    #
    # In anesthetic 2.x the DataFrame uses a MultiIndex of (name, label) tuples.
    # PolyChord parameter columns have integer names (0, 1, ...) when no
    # tidal.paramnames file is present; anesthetic's own columns ('logL',
    # 'logL_birth', 'nlive') have string names.  We identify parameters as
    # columns whose first MultiIndex level is an integer.
    cols = list(ns.columns)
    if cols and isinstance(cols[0], tuple):
        # anesthetic 2.x MultiIndex — strip label level
        param_col_indices = [i for i, c in enumerate(cols) if isinstance(c[0], int)]
    else:
        # anesthetic 1.x or flat index — all columns are parameters
        param_col_indices = list(range(len(cols)))

    ndim = len(param_col_indices)

    # User-supplied names take priority; else fall back to tidal.paramnames;
    # else use p0, p1, ... placeholders.
    if args.param_names:
        user_names = [n.strip() for n in args.param_names.split(",")]
        if len(user_names) != ndim:
            sys.exit(
                f"error: --param-names has {len(user_names)} entries but chain "
                f"has {ndim} parameter columns"
            )
        param_names = user_names
    else:
        from_file = _read_paramnames(chain_root)
        if from_file is not None and len(from_file) == ndim:
            param_names = from_file
        else:
            param_names = [f"p{i}" for i in range(ndim)]
            print(
                f"warning: no tidal.paramnames found; using placeholder names "
                f"{param_names}. Pass --param-names to label them correctly."
            )

    print(f"Parameters ({ndim}): {param_names}")

    # Extract samples — select only the integer-indexed parameter columns
    samples = ns.iloc[:, param_col_indices].to_numpy(dtype=float)
    logl = np.asarray(ns.logL)
    weights = np.asarray(ns.get_weights())

    logz = float(ns.logZ())
    logz_err = float(ns.logZ(100).std())

    # Count dead points: those where logL_birth < logL (live points have logL_birth = -inf)
    n_dead = int(np.isfinite(np.asarray(ns.logL_birth)).sum())

    print(f"Dead points: {n_dead}")
    print(f"Total samples (dead+live): {len(samples)}")
    print(f"log Z = {logz:.3f} ± {logz_err:.3f}")

    # Build InferenceResult and save
    # Import here so this script works standalone without installing tidal as a package
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from tidal.inference._results import InferenceResult

    result = InferenceResult(
        samples=samples,
        log_likelihood=logl,
        log_prior=np.zeros(len(logl)),
        param_names=param_names,
        method="nested",
        log_evidence=logz,
        log_evidence_err=logz_err,
        weights=weights,
        metadata={
            "sampler": "polychord",
            "chain_root": chain_root,
            "_anesthetic_samples": ns,
            "reconstructed": True,
            "n_dead_pts": n_dead,
            "label": args.label,
        },
    )

    result.save(output_path)

    # Patch inference.json to add reconstruction flags (save() doesn't expose metadata keys)
    inf_path = output_path / "inference.json"
    with inf_path.open() as f:
        summary = json.load(f)
    summary["reconstructed"] = True
    summary["n_dead_pts"] = n_dead
    if args.label:
        summary["label"] = args.label
    with inf_path.open("w") as f:
        json.dump(summary, f, indent=2)

    ess = summary.get("effective_sample_size", "?")
    print(f"ESS: {ess}")
    print(f"Written: {output_path}/inference.json")
    print(f"Written: {output_path}/results.csv")
    print()
    print("Now run:")
    print(
        f"  uv run tidal plot {output_path} --type corner --output {output_path}/corner_reconstructed.png"
    )


if __name__ == "__main__":
    main()
