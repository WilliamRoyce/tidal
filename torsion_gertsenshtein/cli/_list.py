"""``tg list`` — List available JSON equation specifications."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

from torsion_gertsenshtein.cli._inspect import discover_parameters

if TYPE_CHECKING:
    from argparse import Namespace


def _find_examples_dir() -> Path:
    """Locate the default examples/data directory."""
    # Try relative to CWD first
    cwd_data = Path.cwd() / "examples" / "data"
    if cwd_data.is_dir():
        return cwd_data

    # Try relative to package
    pkg = Path(__file__).parent.parent
    pkg_data = pkg.parent / "examples" / "data"
    if pkg_data.is_dir():
        return pkg_data

    return cwd_data  # Will fail gracefully in list_command


def list_command(args: Namespace) -> int:
    """Execute the list command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    from torsion_gertsenshtein.symbolic.json_loader import load_equation_system

    scan_dir = Path(args.dir) if args.dir else _find_examples_dir()

    if not scan_dir.is_dir():
        print(f"Error: directory not found: {scan_dir}", file=sys.stderr)
        return 1

    json_files = sorted(scan_dir.glob("*.json"))
    if not json_files:
        print(f"No JSON files found in {scan_dir}")
        return 0

    print(f"Available equation specifications ({scan_dir}):")
    print()

    # Compute column widths for alignment
    entries: list[tuple[str, str, str, str]] = []
    for jf in json_files:
        try:
            spec = load_equation_system(jf)
        except (ValueError, KeyError, OSError) as exc:
            entries.append((jf.name, "?", "?", f"load error: {type(exc).__name__}"))
            continue

        dim_label = f"{spec.spatial_dimension}+1D"
        n_fields = spec.n_components
        field_word = "field" if n_fields == 1 else "fields"
        field_str = f"{n_fields} {field_word}"

        params = discover_parameters(spec)
        if params:
            param_str = "params: " + ", ".join(sorted(params))
        else:
            param_str = "no params needed"

        entries.append((jf.name, dim_label, field_str, param_str))

    # Print aligned
    name_w = max(len(e[0]) for e in entries) + 2
    dim_w = max(len(e[1]) for e in entries) + 2
    field_w = max(len(e[2]) for e in entries) + 2

    for name, dim, fields, params in entries:
        print(f"  {name:<{name_w}} {dim:<{dim_w}} {fields:<{field_w}} {params}")

    print()
    print(f"{len(entries)} specification{'s' if len(entries) != 1 else ''} found.")

    return 0
