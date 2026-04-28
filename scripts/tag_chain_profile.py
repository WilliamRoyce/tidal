#!/usr/bin/env python3
"""Retroactively tag a saved inference chain with a ``stability_profile``.

In-flight chains submitted before the profile registry landed
(commits prior to the #323 resolution) do not carry a
``stability_profile`` field in their ``inference.json`` manifest.
The post-hoc audit needs the field to know which probe was applied;
without it the audit cannot judge whether a sample was
borderline-stable under the historical guard.

This script writes the field in place. Use it once per chain after
``hpc_shuttle.sh pull <jobid>`` lands the directory locally.

Usage
-----

    python scripts/tag_chain_profile.py PATH [--profile NAME] [--force]

* ``PATH``      — directory containing ``inference.json``
                 (e.g. ``hpc_results/ricciem_amplify_28544894``).
* ``--profile`` — profile name. Defaults to ``v1-unit-0.3`` (legacy
                 probe, default for v5/D1/D2/D3).
* ``--force``   — overwrite an existing ``stability_profile`` field.

Examples
--------
    # Tag the in-flight D2.0 chains with the legacy probe.
    python scripts/tag_chain_profile.py hpc_results/bahamonde_amplify_28544894
    python scripts/tag_chain_profile.py hpc_results/bahamonde_suppress_28544902

    # Tag a chain that ran under the new (validated) probe.
    python scripts/tag_chain_profile.py hpc_results/x --profile v2-consistent-0.1

References
----------
- Issue #323 (probe profile registry resolution).
- ``tidal/measurement/_stability_profile.py`` (registered profiles).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("path", type=Path, help="Inference output directory")
    parser.add_argument(
        "--profile",
        default="v1-unit-0.3",
        help="Profile name to tag (default: v1-unit-0.3, the legacy probe)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing stability_profile field",
    )
    args = parser.parse_args(argv)

    # Validate the profile is registered. Importing the registry only after
    # arg parsing keeps `--help` cheap.
    from tidal.measurement._stability_profile import PROFILES

    if args.profile not in PROFILES:
        print(
            f"error: unknown profile '{args.profile}'. Registered: {sorted(PROFILES)}",
            file=sys.stderr,
        )
        return 1

    manifest_path = args.path / "inference.json"
    if not manifest_path.exists():
        print(
            f"error: no inference.json in {args.path}; "
            f"is it an inference output directory?",
            file=sys.stderr,
        )
        return 1

    with manifest_path.open(encoding="utf-8") as f:
        manifest = json.load(f)

    existing = manifest.get("stability_profile")
    if existing and not args.force:
        if existing == args.profile:
            print(f"{manifest_path}: already tagged {existing!r}; nothing to do.")
            return 0
        print(
            f"error: {manifest_path} already tagged {existing!r}; "
            f"refuse to overwrite without --force.",
            file=sys.stderr,
        )
        return 1

    manifest["stability_profile"] = args.profile
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"{manifest_path}: tagged stability_profile={args.profile!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
