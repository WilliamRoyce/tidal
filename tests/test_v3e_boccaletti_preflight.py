"""Tests for the Phase E pre-flight Boccaletti-safety check."""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from typing import TYPE_CHECKING  # noqa: E402

import v3e_boccaletti_preflight as preflight  # noqa: E402

if TYPE_CHECKING:
    import pytest


def _write_env(tmp_path: Path, *, bpeak: float, sigb: float) -> Path:
    env = tmp_path / "_geometry.env"
    env.write_text(
        f"export BPEAK={bpeak}\nexport SIGB={sigb}\n",
    )
    return env


def test_default_geometry_passes() -> None:
    """The canonical Phase E geometry must pass preflight."""
    assert preflight.check(preflight.DEFAULT_ENV) == 0


def test_perturbative_arg_is_safe(tmp_path: Path) -> None:
    env = _write_env(tmp_path, bpeak=0.01, sigb=5.0)
    metrics = preflight.evaluate(preflight._parse_env(env))
    assert 1e-3 <= metrics["arg"] <= 0.3
    assert metrics["dist_to_nearest_higher_zero"] > 0.3
    assert preflight.check(env) == 0


def test_near_higher_zero_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """At Bpeak*sigB ~ 2*pi/sqrt(2*pi) we hit the first higher zero arg = pi."""
    # Solve kappa * Bpeak * sigB * sqrt(2*pi) / 2 = pi  => Bpeak * sigB = 2*pi/sqrt(2*pi).
    target = 2 * math.pi / math.sqrt(2 * math.pi)
    env = _write_env(tmp_path, bpeak=target / 5.0, sigb=5.0)
    assert preflight.check(env) != 0
    out = capsys.readouterr().out
    assert "sin^2 zero" in out


def test_above_perturbative_window_fails(tmp_path: Path) -> None:
    env = _write_env(tmp_path, bpeak=1.0, sigb=5.0)  # arg ~ 6.27 -> well outside
    assert preflight.check(env) != 0


def test_below_perturbative_window_fails(tmp_path: Path) -> None:
    env = _write_env(tmp_path, bpeak=1e-9, sigb=5.0)  # arg ~ 6.27e-9 -> below floor
    assert preflight.check(env) != 0


def test_missing_env_returns_error(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.env"
    rc = preflight.main(["preflight", str(missing)])
    assert rc == 2
