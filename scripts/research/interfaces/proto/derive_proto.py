"""Prototype of the user-facing ``tidalcosmo derive`` command (R-1, #566).

Usage (mirrors ``cobaya-install input.yaml`` then ``cobaya-run input.yaml``)::

    python derive_proto.py runs/vector_one.yaml            # derive if not already stored
    python derive_proto.py runs/vector_one.yaml --test     # check only; never starts a kernel
    python derive_proto.py theories/VectorTheory.yaml --dry-run
    python derive_proto.py theories/bad/x.yaml --validate-only   # stop after Wolfram validation

Design evidence: ``docs/cosmology/r1_planning_record.md`` §3.2-§3.3. The user never types
``wolframscript``; this script launches the FIXED driver as one subprocess, passes the theory
as WXF *data* (kernel-free ``wolframclient`` serializer), and judges success from the
artifact and sentinel lines, never from the exit status alone.
"""

from __future__ import annotations

# cspell:words waveop
import argparse
import json
import os
import platform
import shutil
import subprocess  # noqa: S404 - launching processes is this script's job
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

from wolframclient.deserializers import binary_deserialize
from wolframclient.serializers import export

sys.path.insert(0, str(Path(__file__).resolve().parent))
import theory_store as ts

HERE = Path(__file__).resolve().parent
DRIVER = HERE.parent / "wolfram" / "driver.wls"
REPO_ROOT = HERE.parents[3]
EXPECTED_ENGINE = (
    Path.home()
    / ".local"
    / "wolfram"
    / "engine"
    / "14.3"
    / "Executables"
    / "WolframKernel"
)


def scrub(text: str) -> str:
    """Remove machine-specific paths so transcripts can be committed."""
    return text.replace(str(REPO_ROOT), "<repo>").replace(str(Path.home()), "<home>")


def lane_busy() -> bool:
    """The lane hook cannot see a Python-launched kernel, so guard here (run_tier1_gate.sh:76-83)."""
    for name in ("wolframscript", "WolframKernel"):
        if (
            subprocess.run(
                ["pgrep", "-x", name], capture_output=True, check=False
            ).returncode
            == 0
        ):
            return True
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "input",
        help="a Cobaya run YAML (theory via model: !defaults ...) or a theory YAML",
    )
    parser.add_argument(
        "--store", help="writable spectra store (default: first of the search path)"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="check only: exit 0 if a verified spectrum exists, else 1",
    )
    parser.add_argument(
        "--force", action="store_true", help="derive even if a verified spectrum exists"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="print what would happen; start nothing"
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="run the Wolfram validation only; store nothing",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=0,
        help="seconds; 0 = no cap (derivations can take hours)",
    )
    parser.add_argument("--log-dir", help="also copy the kernel log here")
    args = parser.parse_args(argv)

    theory = ts.load_theory(args.input)
    fp = ts.fingerprint(theory)
    entry, stores = ts.find_entry(fp, args.store)
    print(f"DERIVE_FINGERPRINT={fp}")
    print(f"DERIVE_COUPLINGS={','.join(ts.couplings(theory))}")

    if entry is not None and not args.force and not args.validate_only:
        try:
            ts.verify_entry(entry, theory)
        except ValueError as err:
            print(f"DERIVE_STORE_INVALID={scrub(str(entry))}: {err}")
        else:
            print(f"DERIVE_UP_TO_DATE={scrub(str(entry))}")
            return 0
    if args.test:
        print(
            "DERIVE_TEST=missing: no verified derived spectrum for this theory; searched "
            + ", ".join(scrub(str(s)) for s in stores)
        )
        return 1

    target_store = Path(args.store) if args.store else stores[0]
    final = ts.entry_dir(target_store, fp)
    command = [
        "wolframscript",
        "-file",
        str(DRIVER),
        "<workdir>/theory.wxf",
        "<workdir>",
    ]
    if args.validate_only:
        command.append("--validate-only")
    if args.dry_run:
        print(
            "DERIVE_DRY_RUN would run: "
            + scrub(" ".join(command))
            + f"  -> {scrub(str(final))}"
        )
        return 0

    if not EXPECTED_ENGINE.is_file():
        print(
            "DERIVE_ERROR=no local Wolfram kernel at the pinned engine path; wolframscript would fall back to the cloud (#559)"
        )
        return 1
    if lane_busy():
        print(
            "DERIVE_ERROR=another Wolfram kernel is running; one kernel at a time, machine-wide"
        )
        return 1

    target_store.mkdir(parents=True, exist_ok=True)
    workdir = target_store / f".tmp-{fp[:12]}-{os.getpid()}"
    workdir.mkdir(parents=True)
    theory_wxf = workdir / "theory.wxf"
    theory_wxf.write_bytes(export({**theory, "fingerprint": fp}, target_format="wxf"))

    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    cmd = ["wolframscript", "-file", str(DRIVER), str(theory_wxf), str(workdir)]
    if args.validate_only:
        cmd.append("--validate-only")
    started = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            cwd=workdir,
            env=env,
            capture_output=True,
            text=True,
            timeout=args.timeout or None,
            check=False,
        )
        output, code = proc.stdout + proc.stderr, proc.returncode
    except subprocess.TimeoutExpired as err:
        output = (
            (err.stdout or "") + (err.stderr or "")
            if isinstance(err.stdout, str)
            else ""
        )
        code = "timeout"
    wall = time.monotonic() - started
    (workdir / "derive.log").write_text(output)
    if args.log_dir:
        Path(args.log_dir).mkdir(parents=True, exist_ok=True)
        (Path(args.log_dir) / f"derive-{fp[:12]}.log").write_text(scrub(output))

    sentinels = [
        ln for ln in output.splitlines() if ln.startswith(("STAGE1", "DERIVE"))
    ]
    for line in sentinels:
        print(scrub(line))
    print(f"DERIVE_EXIT_STATUS={code}")
    print(f"DERIVE_KERNEL_WALL_S={wall:.1f}")

    if args.validate_only:
        ok = "STAGE1 PASSED" in output
        shutil.rmtree(workdir)
        return 0 if ok else 1

    spectrum = workdir / "spectrum.wxf"
    # Verdict from the artifact, never from the exit status (Wolfram 14.3 shutdown segfaults; PSALTer's bare Quit[]).
    if "STAGE1 PASSED" not in output or not spectrum.is_file():
        print("DERIVE_ERROR=no verified spectrum was produced; nothing stored")
        shutil.rmtree(workdir)
        return 1
    data = binary_deserialize(spectrum.read_bytes())
    if data.get("fingerprint") != fp:
        print(
            "DERIVE_ERROR=spectrum.wxf carries a different fingerprint; nothing stored"
        )
        shutil.rmtree(workdir)
        return 1

    manifest = {
        "schema": ts.SCHEMA,
        "fingerprint": fp,
        "theory": theory,
        "couplings": ts.couplings(theory),
        "conventions": ts.CONVENTIONS,
        "spectrum_sha256": ts.sha256_file(spectrum),
        "generator": {
            "psalter_version": data.get("psalter_version"),
            "wolfram_version": str(data.get("wolfram_version")),
            "driver": "scripts/research/interfaces/wolfram/driver.wls",
        },
        "waveop_sha256": next(
            (
                ln.split("=", 1)[1]
                for ln in sentinels
                if ln.startswith("STAGE1_WAVEOP_SHA256=")
            ),
            None,
        ),
        "particle_spectrum_wall_s": data.get("particle_spectrum_wall_s"),
        "kernel_wall_s": round(wall, 1),
        "host": platform.node(),
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    (workdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True)
    )
    theory_wxf.unlink()
    for extra in (
        workdir.iterdir()
    ):  # PSALTer's own .mx/.pdf go beside the log, not into the entry
        if extra.name not in {"spectrum.wxf", "manifest.json", "derive.log"}:
            (workdir / "psalter_outputs").mkdir(exist_ok=True)
            extra.rename(workdir / "psalter_outputs" / extra.name)
    final.parent.mkdir(parents=True, exist_ok=True)
    if final.exists():
        shutil.rmtree(final)
    workdir.rename(final)  # atomic on one filesystem: readers never see a partial entry
    ts.verify_entry(final, theory)
    print(f"DERIVE_STORED={scrub(str(final))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
