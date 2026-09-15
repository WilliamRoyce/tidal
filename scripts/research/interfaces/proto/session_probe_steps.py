"""R-1 prototype: isolate the ``wolframclient`` session measurements (follow-up to session_probe.py).

The first probe hung inside the PSALTer derivation for its full 600 s watchdog, which left
the kernel busy, so its timeout and ``Quit[]`` measurements were confounded. This probe
separates them, each in a FRESH session whose kernels are killed before the next starts:

* ``bare``: no PSALTer. Is ``evaluate_wrap(..., timeout=2)`` honored on an 8 s evaluation
  (wolframclient #38)? What does a bare ``Quit[]`` do to the session?
* ``steps``: the Stage-1 route one step at a time, each with its own watchdog — ``Needs``,
  ``Get`` of the package, ``RunStage1`` with ``ValidateOnly``, then ``DefField`` alone — to
  locate where a session run stops making progress.

Usage::

    python session_probe_steps.py bare  <workdir>
    python session_probe_steps.py steps <theory.yaml> <workdir>
"""

from __future__ import annotations

import os
import signal
import subprocess  # noqa: S404 - process inspection is this script's job
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
    print(f"STEPS_{key}={text}", flush=True)


def watchdog(fn, seconds: float) -> tuple[str, object, float]:
    box: dict[str, object] = {}

    def target() -> None:
        try:
            box["value"], box["outcome"] = fn(), "ok"
        except BaseException as err:
            box["value"], box["outcome"] = f"{type(err).__name__}: {err}", "raised"

    start = time.monotonic()
    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    thread.join(seconds)
    elapsed = time.monotonic() - start
    return (
        ("hung", None, elapsed)
        if thread.is_alive()
        else (box["outcome"], box.get("value"), elapsed)
    )


def kill_kernels() -> int:
    out = subprocess.run(
        ["pgrep", "-x", "WolframKernel"], capture_output=True, text=True, check=False
    ).stdout
    pids = [int(p) for p in out.split()]
    for pid in pids:
        os.kill(pid, signal.SIGKILL)
    time.sleep(1)
    return len(pids)


def fresh_session(workdir: Path) -> WolframLanguageSession:
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    os.chdir(workdir)
    session = WolframLanguageSession(kernel=str(KERNEL))
    outcome, value, elapsed = watchdog(session.start, 120)
    say("START", f"{outcome} {elapsed:.1f}s {value or ''}")
    return session


def bare(workdir: Path) -> None:
    session = fresh_session(workdir)
    outcome, value, elapsed = watchdog(
        lambda: session.evaluate_wrap("Pause[8]; 1", timeout=2), 30
    )
    say(
        "BARE_TIMEOUT",
        f"timeout=2 on Pause[8]: {outcome} after {elapsed:.1f}s ({value})",
    )
    outcome, value, elapsed = watchdog(lambda: session.evaluate("1+1"), 30)
    say("BARE_AFTER_TIMEOUT_1PLUS1", f"{outcome} after {elapsed:.1f}s ({value})")
    outcome, value, elapsed = watchdog(lambda: session.evaluate("Quit[]"), 30)
    say("BARE_QUIT", f"{outcome} after {elapsed:.1f}s ({value})")
    outcome, value, elapsed = watchdog(lambda: session.evaluate("1+1"), 20)
    say("BARE_AFTER_QUIT_1PLUS1", f"{outcome} after {elapsed:.1f}s ({value})")
    outcome, value, elapsed = watchdog(session.terminate, 20)
    say("BARE_TERMINATE", f"{outcome} after {elapsed:.1f}s ({value})")
    say("BARE_KERNELS_KILLED_AFTERWARDS", kill_kernels())


def steps(theory_path: Path, workdir: Path) -> None:
    theory = ts.load_theory(theory_path)
    fp = ts.fingerprint(theory)
    (workdir / "theory.wxf").write_bytes(
        export({**theory, "fingerprint": fp}, target_format="wxf")
    )
    session = fresh_session(workdir)
    plan = [
        ("SETDIR", f'SetDirectory["{workdir}"]', 20),
        ("NEEDS_PSALTER", 'Needs["xAct`PSALTer`"]; "loaded"', 90),
        ("GET_PACKAGE", f'Get["{PACKAGE}"]; "got"', 60),
        (
            "VALIDATE_ONLY",
            f'R1Proto`RunStage1[Import["{workdir / "theory.wxf"}", "WXF"], "{workdir}", "ValidateOnly" -> True]',
            60,
        ),
        (
            "DEFFIELD_ALONE",
            'DefField[ProbeField[-a], PrintAs -> "B", PrintSourceAs -> "k"]; "defined"',
            120,
        ),
    ]
    for label, code, limit in plan:
        outcome, value, elapsed = watchdog(
            lambda c=code: session.evaluate_wrap(c), limit
        )
        detail = (
            f"success={value.success} result={value.result}"
            if outcome == "ok"
            else value
        )
        say(label, f"{outcome} after {elapsed:.1f}s ({detail})")
        if outcome != "ok":
            break
    say("FILES_IN_WORKDIR", sorted(p.name for p in workdir.iterdir()))
    say("KERNELS_KILLED_AFTERWARDS", kill_kernels())


def main() -> int:
    out = subprocess.run(
        ["pgrep", "-x", "WolframKernel"], capture_output=True, check=False
    ).returncode
    if (
        out == 0
        or subprocess.run(
            ["pgrep", "-x", "wolframscript"], capture_output=True, check=False
        ).returncode
        == 0
    ):
        say("ERROR", "another Wolfram kernel is running; one kernel at a time")
        return 1
    mode = sys.argv[1]
    workdir = Path(sys.argv[-1]).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    if mode == "bare":
        bare(workdir)
    else:
        steps(Path(sys.argv[2]).resolve(), workdir)
    return 0


if __name__ == "__main__":
    code = main()
    os._exit(code)
