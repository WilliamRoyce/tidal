"""R-1 prototype: measure the ``wolframclient`` session route instead of arguing about it.

Runs the SAME committed package (``wolfram/Stage1Proto.wl``) inside a
``WolframLanguageSession`` and records three things the plan's rejection of sessions rests on
(``docs/cosmology/r1_planning_record.md`` §2.4):

1. whether the package derives the Vector theory through a session at all (spectrum.wxf is
   written for the wave-operator comparison);
2. whether a Python-side ``timeout=`` is honored (wolframclient #38 says it is ignored);
3. what a bare ``Quit[]`` — PSALTer has one at ``ConjectureInverse.m:37-39`` — does to the
   session: exception, hang, or silent loss.

Every step has a watchdog, and any kernel this script started is killed before it exits, so
the machine-wide single-kernel lane is never left held. Usage::

    python session_probe.py <theory.yaml> <workdir>
"""

from __future__ import annotations

import os
import signal
import subprocess  # noqa: S404 - launching processes is this script's job
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import theory_store as ts
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.serializers import export

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent / "wolfram" / "Stage1Proto.wl"
KERNEL = (
    Path.home()
    / ".local"
    / "wolfram"
    / "engine"
    / "14.3"
    / "Executables"
    / "WolframKernel"
)
REPO_ROOT = HERE.parents[3]


def say(key: str, value: object) -> None:
    text = (
        str(value).replace(str(REPO_ROOT), "<repo>").replace(str(Path.home()), "<home>")
    )
    print(f"SESSION_{key}={text}", flush=True)


def with_watchdog(fn, seconds: float):
    """Run fn in a thread; return (outcome, value, elapsed) with outcome ok|raised|hung."""
    box: dict[str, object] = {}

    def target() -> None:
        try:
            box["value"] = fn()
            box["outcome"] = "ok"
        except BaseException as err:
            box["value"] = f"{type(err).__name__}: {err}"
            box["outcome"] = "raised"

    start = time.monotonic()
    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    thread.join(seconds)
    elapsed = time.monotonic() - start
    if thread.is_alive():
        return "hung", None, elapsed
    return box["outcome"], box.get("value"), elapsed


def kernel_pids() -> set[int]:
    out = subprocess.run(
        ["pgrep", "-x", "WolframKernel"], capture_output=True, text=True, check=False
    ).stdout
    return {int(p) for p in out.split()}


def main() -> int:
    theory_path, workdir = Path(sys.argv[1]), Path(sys.argv[2])
    if (
        kernel_pids()
        or subprocess.run(
            ["pgrep", "-x", "wolframscript"], capture_output=True, check=False
        ).returncode
        == 0
    ):
        say("ERROR", "another Wolfram kernel is running; one kernel at a time")
        return 1
    workdir.mkdir(parents=True, exist_ok=True)
    theory = ts.load_theory(theory_path)
    fp = ts.fingerprint(theory)
    (workdir / "theory.wxf").write_bytes(
        export({**theory, "fingerprint": fp}, target_format="wxf")
    )

    os.environ["QT_QPA_PLATFORM"] = (
        "offscreen"  # inherited by the kernel the session spawns
    )
    os.chdir(workdir)  # PSALTer freezes $WorkingDirectory at load
    session = WolframLanguageSession(kernel=str(KERNEL))

    outcome, value, elapsed = with_watchdog(session.start, 120)
    say("START", f"{outcome} {elapsed:.1f}s {value if outcome != 'ok' else ''}")
    if outcome != "ok":
        return 1

    # 1. The same package through the session.
    def derive():
        session.evaluate(f'SetDirectory["{workdir}"]')
        session.evaluate('Needs["xAct`PSALTer`"]')
        session.evaluate(f'Get["{PACKAGE}"]')
        return session.evaluate_wrap(
            f'R1Proto`RunStage1[Import["{workdir / "theory.wxf"}", "WXF"], "{workdir}"]'
        )

    outcome, result, elapsed = with_watchdog(derive, 600)
    say("DERIVE_OUTCOME", f"{outcome} after {elapsed:.1f}s")
    if outcome == "ok":
        say("DERIVE_SUCCESS", result.success)
        say("DERIVE_RESULT", result.result)
        say("DERIVE_MESSAGE_COUNT", len(result.messages or []))
        say(
            "DERIVE_MESSAGE_NAMES",
            sorted({str(n) for n in (result.messages_name or [])})[:12],
        )
    say("DERIVE_SPECTRUM_WRITTEN", (workdir / "spectrum.wxf").is_file())

    # 2. Is a Python-side timeout honored?
    outcome, value, elapsed = with_watchdog(
        lambda: session.evaluate_wrap("Pause[8]; 1", timeout=2), 60
    )
    say(
        "TIMEOUT_PROBE",
        f"asked timeout=2s on an 8s evaluation: {outcome} after {elapsed:.1f}s ({value})",
    )

    # 3. A bare Quit[] inside the session (PSALTer's ConjectureInverse.m:37-39 does this).
    before = kernel_pids()
    outcome, value, elapsed = with_watchdog(lambda: session.evaluate("Quit[]"), 60)
    say("QUIT_PROBE", f"{outcome} after {elapsed:.1f}s ({value})")
    outcome, value, elapsed = with_watchdog(lambda: session.evaluate("1+1"), 30)
    say("AFTER_QUIT_EVALUATE", f"{outcome} after {elapsed:.1f}s ({value})")
    outcome, value, elapsed = with_watchdog(session.terminate, 30)
    say("TERMINATE", f"{outcome} after {elapsed:.1f}s ({value})")

    leftover = kernel_pids()
    say("KERNELS_BEFORE_QUIT", len(before))
    say("KERNELS_LEFT_AFTER_TERMINATE", len(leftover))
    for pid in leftover:  # never leave the lane held
        os.kill(pid, signal.SIGKILL)
    time.sleep(1)
    say("KERNELS_AFTER_CLEANUP", len(kernel_pids()))
    return 0


if __name__ == "__main__":
    code = main()
    os._exit(code)  # a hung watchdog thread must not keep the process alive
