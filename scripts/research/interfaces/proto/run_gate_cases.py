"""R-1 prototype: exercise the Cobaya gate (``spectator_gate.py``) — refusing, then passing.

No Wolfram. Uses the derived spectrum already stored by ``derive_proto.py``. Each case runs
Cobaya in its own subprocess from real YAML files (so ``!defaults`` includes are exercised
exactly as ``cobaya-run`` does), and prints a scrubbed transcript. Usage::

    python run_gate_cases.py <store-with-derived-vector-spectrum> <scratch-dir>
"""

from __future__ import annotations

# cspell:words spawnl spawnv execv
import ast
import json
import shutil
import subprocess  # noqa: S404 - launching processes is this script's job
import sys
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO_ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
import theory_store as ts

PARAMS_OK = """
params:
  Theta1: {prior: {min: 0, max: 2}, ref: 1.0, proposal: 0.05, latex: \\theta_1}
  Theta2: {prior: {min: 0, max: 2}, latex: \\theta_2}
  Theta3: 0.5
"""


def scrub(text: str, *extra: Path) -> str:
    for p in (*extra, REPO_ROOT, Path.home()):
        text = text.replace(
            str(p),
            f"<{p.name if p in extra else ('repo' if p == REPO_ROOT else 'home')}>",
        )
    return text


def run_file(
    yaml_path: Path, env_store: str, resume_output: str | None = None
) -> tuple[int, str]:
    """get_model + logposterior (or cobaya.run with output) in a fresh interpreter."""
    code = textwrap.dedent(f"""
        import sys
        from cobaya.model import get_model
        from cobaya.run import run
        if {resume_output!r}:
            info = {{}}
            upd, sampler = run({str(yaml_path)!r}, output={resume_output!r}, resume=True)
            print("RUN_OK")
        else:
            model = get_model({str(yaml_path)!r})
            point = model.prior.reference()  # whatever is sampled in this run file
            print("SAMPLED", list(model.parameterization.sampled_params()))
            print("LOGPOST", model.logposterior(point).logpost)
            print("VERSION", model.theory["spectator_gate.SpectatorGate"].get_version()[:21])
        """)
    env = {
        "TIDALCOSMO_SPECTRA_PATH": env_store,
        "PATH": "/usr/bin:/bin",
        "HOME": str(Path.home()),
    }
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
        check=False,
        cwd=yaml_path.parent,
    )
    return proc.returncode, proc.stdout + proc.stderr


def key_lines(out: str) -> str:
    keep = [
        ln
        for ln in out.splitlines()
        if any(
            k in ln
            for k in (
                "SAMPLED",
                "LOGPOST",
                "VERSION",
                "RUN_OK",
                "*ERROR*",
                "Error",
                "error",
                "Requirement",
                "Could not find",
                "not compatible",
            )
        )
    ]
    return "\n".join(dict.fromkeys(keep)) or out[-600:]


ONE = "likelihood:\n  one:\n"
NON_ABSORBING = "likelihood:\n  flat:\n    external: 'lambda: 0.0'\n"
"""A likelihood that claims no parameters: unlike ``one`` (an AbsorbUnusedParamsLikelihood), it
does not silently swallow a misspelled parameter, so Cobaya's own typo error can fire."""


def write_run(
    scratch: Path,
    name: str,
    theory_file: Path,
    params: str = PARAMS_OK,
    extra_theory: str = "",
    likelihood: str = ONE,
) -> Path:
    path = scratch / f"{name}.yaml"
    path.write_text(
        textwrap.dedent(f"""\
        theory:
          spectator_gate.SpectatorGate:
            python_path: {HERE}
            model: !defaults {theory_file}
        {extra_theory}""")
        + params
        + likelihood
        + "sampler:\n  evaluate:\n"
    )
    return path


def static_no_kernel_check() -> str:
    """The Cobaya-facing modules must contain no way to start a kernel."""
    # The kernel-free WXF (de)serializers are allowed: the production Theory reads spectrum.wxf with
    # wolframclient.deserializers. What is forbidden is anything that can start a process or a session.
    forbidden_modules = {
        "subprocess",
        "multiprocessing",
        "pexpect",
        "wolframclient.evaluation",
    }
    forbidden_calls = {
        "system",
        "popen",
        "spawnl",
        "spawnv",
        "execv",
        "execl",
        "startfile",
    }
    findings = []
    for mod in ("spectator_gate.py", "theory_store.py"):
        tree = ast.parse((HERE / mod).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                findings += [
                    f"{mod}: import {a.name}"
                    for a in node.names
                    if a.name in forbidden_modules
                    or a.name.startswith("wolframclient.evaluation")
                ]
            elif isinstance(node, ast.ImportFrom) and node.module:
                if node.module in forbidden_modules or node.module.startswith(
                    "wolframclient.evaluation"
                ):
                    findings.append(f"{mod}: from {node.module} import ...")
            elif isinstance(node, ast.Attribute) and node.attr in forbidden_calls:
                findings.append(f"{mod}: .{node.attr}")
    return (
        "PASS (no kernel launcher reachable)"
        if not findings
        else "FAIL " + "; ".join(findings)
    )


def main() -> int:
    store, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    theory_file = ROOT / "theories" / "VectorTheory.yaml"
    theory = ts.load_theory(theory_file)
    fp = ts.fingerprint(theory)
    entry = ts.entry_dir(store, fp)
    assert (entry / "manifest.json").is_file(), (
        f"derive the Vector theory into {store} first"
    )
    empty = scratch / "empty_store"
    empty.mkdir()

    def case(
        label: str, yaml_path: Path, env_store: Path | str, expect_ok: bool, **kw: str
    ) -> None:
        code, out = run_file(yaml_path, str(env_store), **kw)
        ok = (code == 0) == expect_ok
        print(
            f"\n=== CASE {label}: exit={code} expected={'pass' if expect_ok else 'refuse'} -> {'AS EXPECTED' if ok else 'UNEXPECTED'}"
        )
        print(scrub(key_lines(out), store, scratch))

    run_ok = write_run(scratch, "run_ok", theory_file)

    # --- the refusals -----------------------------------------------------------------
    case(
        "1 no derived spectrum anywhere on the search path",
        run_ok,
        empty,
        expect_ok=False,
    )

    edited = scratch / "VectorTheory_edited.yaml"
    edited.write_text(
        theory_file.read_text().replace("Theta3: '-(1/2)*", "Theta3: '-(1/3)*")
    )
    case(
        "2 operator edited after derivation (fingerprint changes)",
        write_run(scratch, "run_edited", edited),
        store,
        expect_ok=False,
    )

    corrupt = scratch / "corrupt_store"
    shutil.copytree(
        store, corrupt, ignore=shutil.ignore_patterns("psalter_outputs", "derive.log")
    )
    corrupt_wxf = ts.entry_dir(corrupt, fp) / "spectrum.wxf"
    corrupt_wxf.write_bytes(corrupt_wxf.read_bytes()[:-8])
    case("3 spectrum.wxf truncated (partial copy)", run_ok, corrupt, expect_ok=False)

    wrongconv = scratch / "wrongconv_store"
    shutil.copytree(
        store, wrongconv, ignore=shutil.ignore_patterns("psalter_outputs", "derive.log")
    )
    manifest_file = ts.entry_dir(wrongconv, fp) / "manifest.json"
    manifest = json.loads(manifest_file.read_text())
    manifest["conventions"] = {"signature": [-1, 1, 1, 1], "epsilon0123": -1}
    manifest_file.write_text(json.dumps(manifest))
    case(
        "4 manifest declares legacy conventions (-,+,+,+), eps=-1",
        run_ok,
        wrongconv,
        expect_ok=False,
    )

    missing = write_run(
        scratch,
        "run_missing_param",
        theory_file,
        params="\nparams:\n  Theta1: {prior: {min: 0, max: 2}}\n  Theta3: 0.5\n",
    )
    case("5 coupling Theta2 missing from params", missing, store, expect_ok=False)

    extra = write_run(
        scratch,
        "run_extra_param",
        theory_file,
        params=PARAMS_OK + "  Theta4: {prior: {min: 0, max: 1}}\n",
        likelihood=NON_ABSORBING,
    )
    case(
        "6 misspelled/extra parameter Theta4 (non-absorbing likelihood)",
        extra,
        store,
        expect_ok=False,
    )
    extra_one = write_run(
        scratch,
        "run_extra_param_one",
        theory_file,
        params=PARAMS_OK + "  Theta4: {prior: {min: 0, max: 1}}\n",
    )
    case(
        "6b same typo, but with the unit likelihood `one`, which absorbs unused parameters (caveat)",
        extra_one,
        store,
        expect_ok=True,
    )

    # --- the passes ---------------------------------------------------------------------
    case("7 everything matches", run_ok, store, expect_ok=True)
    reprior = write_run(
        scratch,
        "run_new_priors",
        theory_file,
        params="\nparams:\n  Theta1: {prior: {min: -3, max: 3}}\n  Theta2: 0.3\n  Theta3: {prior: {min: 0, max: 1}}\n",
    )
    case(
        "8 different priors and a different fixed coupling (no re-derivation)",
        reprior,
        store,
        expect_ok=True,
    )
    case(
        "9 found in the second store of the search path",
        run_ok,
        f"{empty}:{store}",
        expect_ok=True,
    )

    # --- resume provenance ----------------------------------------------------------------
    chains = scratch / "chains" / "vector"
    code, out = run_file(run_ok, str(store), resume_output=str(chains))
    print(f"\n=== RESUME setup: first run with output exit={code}")
    print(scrub(key_lines(out), store, scratch))
    upd = next(chains.parent.glob("*.updated.yaml"), None)
    if upd:
        print(
            "updated.yaml records:",
            [
                ln.strip()
                for ln in upd.read_text().splitlines()
                if "version:" in ln and "spectrum-" in ln
            ],
        )
    other = scratch / "OtherTheory.yaml"
    other.write_text(
        theory_file.read_text().replace("Theta3: '-(1/2)*", "Theta3: '-(1/4)*")
    )
    other_store = scratch / "other_store"
    shutil.copytree(
        store,
        other_store,
        ignore=shutil.ignore_patterns("psalter_outputs", "derive.log"),
    )
    other_theory = ts.load_theory(other)
    ofp = ts.fingerprint(other_theory)
    src = ts.entry_dir(other_store, fp)
    dst = ts.entry_dir(other_store, ofp)
    dst.parent.mkdir(parents=True, exist_ok=True)
    src.rename(dst)
    om = json.loads((dst / "manifest.json").read_text())
    om.update(fingerprint=ofp, theory=other_theory)
    (dst / "manifest.json").write_text(json.dumps(om))
    case(
        "10 resume the same chain with a different theory/spectrum",
        write_run(scratch, "run_other", other),
        other_store,
        expect_ok=False,
        resume_output=str(chains),
    )
    moved = scratch / "moved_store"
    shutil.copytree(store, moved)
    case(
        "11 resume the same chain with the SAME spectrum found in a moved store",
        run_ok,
        moved,
        expect_ok=True,
        resume_output=str(chains),
    )

    print("\n=== STATIC no-kernel check:", static_no_kernel_check())
    return 0


if __name__ == "__main__":
    sys.exit(main())
