"""The shared ``run_status`` vocabulary (GH #480).

``run_status`` is the column a reader of a sweep ``results.csv`` or a
chain CSV uses to decide whether a row is usable and what went wrong when
it is not.  It had no single definition until GH #480: ``tidal sweep``
and ``tidal sample`` each grew their own vocabulary, sharing only
``success``, and three documents described three mutually inconsistent
taxonomies -- two of them naming tags no code ever emitted.

The predecessor of this file asserted that each of nine hardcoded tag
strings ``isinstance(tag, str)``.  That test could not fail, and it
listed ``simulation_diverged`` -- a tag nothing emitted -- which is
precisely how the drift survived unnoticed.  The coverage test below is
its replacement: it reads the strings out of the source and checks them
against the enumeration, so a new tag that is not declared fails the
suite.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tidal.measurement._run_stages import RunStatus

TIDAL_ROOT = Path(__file__).resolve().parents[1] / "tidal"


_KEY = "run_status"


def _is_key(node: ast.expr | None) -> bool:
    return isinstance(node, ast.Constant) and node.value == _KEY


def _const_str(node: ast.expr | None) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _from_dict(node: ast.AST) -> list[str]:
    """Catch ``{"run_status": "success"}``."""
    if not isinstance(node, ast.Dict):
        return []
    return [
        v
        for key, val in zip(node.keys, node.values, strict=True)
        if _is_key(key) and (v := _const_str(val)) is not None
    ]


def _from_subscript_assign(node: ast.AST) -> list[str]:
    """Catch ``metrics["run_status"] = "success"``."""
    if not isinstance(node, ast.Assign):
        return []
    value = _const_str(node.value)
    if value is None:
        return []
    return [
        value
        for tgt in node.targets
        if isinstance(tgt, ast.Subscript) and _is_key(tgt.slice)
    ]


def _from_setdefault(node: ast.AST) -> list[str]:
    """Catch ``metrics.setdefault("run_status", "success")``."""
    if not (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "setdefault"
        and len(node.args) == 2
        and _is_key(node.args[0])
    ):
        return []
    value = _const_str(node.args[1])
    return [value] if value is not None else []


def _literals_in_source(source: str, label: str) -> dict[str, set[str]]:
    """String literals assigned to a ``run_status`` key in one source.

    Walks the AST rather than grepping, so every assignment form is
    caught while comments and docstrings mentioning a tag are not.
    """
    found: dict[str, set[str]] = {}
    tree = ast.parse(source, filename=label)
    for node in ast.walk(tree):
        for extract in (_from_dict, _from_subscript_assign, _from_setdefault):
            for value in extract(node):
                found.setdefault(value, set()).add(label)
    return found


def _assigned_run_status_literals() -> dict[str, set[str]]:
    """Every string literal assigned to ``run_status`` across ``tidal/``."""
    found: dict[str, set[str]] = {}
    for path in TIDAL_ROOT.rglob("*.py"):
        label = str(path.relative_to(TIDAL_ROOT))
        for tag, labels in _literals_in_source(
            path.read_text(encoding="utf-8"), label
        ).items():
            found.setdefault(tag, set()).update(labels)
    return found


# Every assignment form the walker must recognize, with a tag that is NOT a
# RunStatus member so a passing coverage test cannot be confused with a
# broken walker.
_WALKER_FIXTURE = '''
def a():
    return {"run_status": "fixture_dict"}
def b(metrics):
    metrics["run_status"] = "fixture_subscript"
def c(metrics):
    metrics.setdefault("run_status", "fixture_setdefault")
def d():
    # "run_status": "fixture_comment" — must NOT be picked up
    """Nor "run_status": "fixture_docstring"."""
'''


class TestRunStatusCoverage:
    """Adding a status without declaring it must fail the suite."""

    def test_the_walker_actually_detects_literals(self) -> None:
        """Prove the detector works before trusting what it does not find.

        Without this, ``test_no_raw_literals`` below would pass just as
        happily against a broken walker — which is how the predecessor of
        this file ended up asserting nothing at all.
        """
        found = _literals_in_source(_WALKER_FIXTURE, "fixture")
        assert set(found) == {
            "fixture_dict",
            "fixture_subscript",
            "fixture_setdefault",
        }, f"walker missed an assignment form: {sorted(found)}"

    def test_no_raw_literals(self) -> None:
        """Statuses must go through :class:`RunStatus`, never a bare string.

        This is what keeps the enumeration authoritative: a raw literal is
        a status that exists in the data but not in the vocabulary, which
        is how sweep and inference drifted apart in the first place.
        """
        literals = _assigned_run_status_literals()
        assert not literals, (
            "raw run_status string literals under tidal/ — use RunStatus.* "
            f"(tidal/measurement/_run_stages.py): "
            f"{ {k: sorted(v) for k, v in literals.items()} }"
        )

    def test_any_literal_that_survives_is_declared(self) -> None:
        """Backstop for a partial regression: a literal is at least known."""
        literals = _assigned_run_status_literals()
        declared = {str(s) for s in RunStatus}
        undeclared = {
            tag: sorted(paths) for tag, paths in literals.items() if tag not in declared
        }
        assert not undeclared, (
            "run_status values emitted but absent from RunStatus "
            f"(tidal/measurement/_run_stages.py): {undeclared}"
        )

    def test_historical_values_are_not_emitted(self) -> None:
        """``diverged`` / ``tachyonic_gated`` are archive-only.

        They stay in the enumeration because recorded CSVs contain them
        and a reader needs to know what they meant -- but nothing should
        write them again.
        """
        literals = set(_assigned_run_status_literals())
        historical = {str(s) for s in RunStatus} - {str(s) for s in RunStatus.live()}
        assert historical == {"diverged", "tachyonic_gated"}
        assert not (literals & historical), (
            f"historical run_status value re-emitted: {literals & historical}"
        )


class TestRunStatusSerialization:
    """The property that lets every existing consumer stay untouched."""

    @pytest.mark.parametrize("status", list(RunStatus))
    def test_compares_equal_to_its_plain_string(self, status: RunStatus) -> None:
        assert status == str(status)

    def test_survives_a_csv_round_trip(self, tmp_path: Path) -> None:
        """Consumers filter rows with ``row["run_status"] == "success"``.

        ``_sweep_panels``, ``examples/*/analyze_sweep.py`` and
        ``scripts/pull_and_plot.sh`` all do this against values that have
        been through a CSV.  A StrEnum writes as its bare value, so they
        need no change -- this pins that.
        """
        import csv

        path = tmp_path / "results.csv"
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["run_status", "P_max"])
            w.writeheader()
            w.writerow({"run_status": RunStatus.SUCCESS, "P_max": 0.5})
            w.writerow({"run_status": RunStatus.SIMULATION_DIVERGED, "P_max": ""})

        with path.open() as f:
            rows = list(csv.DictReader(f))

        assert rows[0]["run_status"] == "success"
        assert rows[1]["run_status"] == "simulation_diverged"
        assert [r for r in rows if r["run_status"] == "success"] == [rows[0]]

    def test_json_serializes_as_a_bare_string(self) -> None:
        import json

        assert json.dumps({"run_status": RunStatus.KINETIC_ERROR}) == (
            '{"run_status": "kinetic_error"}'
        )


class TestExceptionClassification:
    """Distinct causes must not collapse into one physics-sounding tag."""

    def test_divergence_is_its_own_status(self) -> None:
        from tidal.solver import SimulationDivergedError

        assert (
            RunStatus.from_exception(SimulationDivergedError("fields blew up"))
            is RunStatus.SIMULATION_DIVERGED
        )

    def test_kinetic_evaluation_error_is_not_a_physics_verdict(self) -> None:
        """GH #447/#480: a missing ``--param`` is a configuration error.

        It used to be recorded as ``diverged`` by the sweep row path --
        the same mislabeling GH #447 prevented at the probe by making this
        a ``RuntimeError`` rather than a ``ValueError``, reintroduced one
        layer up.
        """
        from tidal.solver import KineticEvaluationError

        status = RunStatus.from_exception(KineticEvaluationError("xi unbound"))
        assert status is RunStatus.KINETIC_ERROR
        assert status is not RunStatus.SIMULATION_DIVERGED

    @pytest.mark.parametrize(
        "exc",
        [
            ValueError("bad"),
            TypeError("bad"),
            KeyError("bad"),
            OSError("bad"),
            RuntimeError("?"),
        ],
    )
    def test_generic_types_are_not_attributed_to_a_stage(self, exc: Exception) -> None:
        """An exception TYPE does not identify the stage that raised it.

        The first cut of this helper mapped ValueError/TypeError/KeyError/
        OSError to MEASUREMENT_ERROR. That over-claims: a bare OSError is a
        missing spec file as readily as an unreadable output, and a
        ValueError is a bad parameter as readily as a failed measurement.
        Concretely it mislabelled a missing spec path on the inference
        path, which the GH #480 equal-output proof caught.

        MEASUREMENT_ERROR is now set by ``run_point`` at the point where it
        is actually measuring — attribution by position, not by type.
        """
        assert RunStatus.from_exception(exc) is RunStatus.EXCEPTION

    def test_measurement_stage_failures_are_attributed_by_position(self) -> None:
        """``run_point`` tags a failure in the measure stage, wherever it came from."""
        from tidal.measurement._run_stages import (
            _MeasurementStageError,  # pyright: ignore[reportPrivateUsage]
        )

        cause = KeyError("P_max")
        wrapper = _MeasurementStageError(cause)
        assert wrapper.cause is cause
        # The type alone would have said EXCEPTION; position says otherwise.
        assert RunStatus.from_exception(cause) is RunStatus.EXCEPTION


class TestRunPointTotality:
    """``run_point`` returns an outcome for every path, and never raises.

    That totality is what lets both callers be straight-line mappings
    instead of each re-implementing exception handling — which is how they
    drifted apart in the first place (GH #454).
    """

    @staticmethod
    def _ctx(tmp_path: Path) -> object:
        from argparse import Namespace

        from tidal.measurement._run_stages import PointContext

        return PointContext(
            spec_path=Path("/nonexistent/spec.json"),
            base_args=Namespace(param=[], grid_shape=None, bounds=None),
            param_overrides={},
            measurements={"summary"},
            source=None,
            target=None,
            threshold=0.99,
            run_dir=tmp_path,
        )

    @pytest.mark.parametrize(
        ("exc_factory", "expected"),
        [
            ("SimulationDivergedError", RunStatus.SIMULATION_DIVERGED),
            ("KineticEvaluationError", RunStatus.KINETIC_ERROR),
            ("ValueError", RunStatus.EXCEPTION),
            ("RuntimeError", RunStatus.EXCEPTION),
        ],
    )
    def test_stage_failures_become_outcomes(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        exc_factory: str,
        expected: RunStatus,
    ) -> None:
        import tidal.measurement._run_stages as stages
        from tidal.solver import KineticEvaluationError, SimulationDivergedError

        kinds: dict[str, type[BaseException]] = {
            "SimulationDivergedError": SimulationDivergedError,
            "KineticEvaluationError": KineticEvaluationError,
            "ValueError": ValueError,
            "RuntimeError": RuntimeError,
        }
        exc_type = kinds[exc_factory]

        def _boom(*_a: object, **_k: object) -> object:
            msg = "stage failed"
            raise exc_type(msg)

        def _no_probe(*_a: object, **_k: object) -> tuple[None, dict[str, object]]:
            return None, {}

        monkeypatch.setattr(stages, "probe_for_run", _no_probe)
        monkeypatch.setattr(stages, "simulate_run", _boom)

        outcome = stages.run_point(self._ctx(tmp_path), backend="disk")

        assert outcome.status is expected
        assert not outcome.ok
        assert isinstance(outcome.exception, exc_type)

    def test_measurement_failure_is_attributed_to_the_measure_stage(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        """A ValueError from measure is MEASUREMENT_ERROR; from simulate it is not.

        Position, not type — the distinction the GH #480 equal-output proof
        forced, after type-sniffing mislabelled a missing spec path.
        """
        import tidal.measurement._run_stages as stages

        def _no_probe(*_a: object, **_k: object) -> tuple[None, dict[str, object]]:
            return None, {}

        def _ok_sim(*_a: object, **_k: object) -> tuple[int, float, None]:
            return 0, 1.0, None

        monkeypatch.setattr(stages, "probe_for_run", _no_probe)
        monkeypatch.setattr(stages, "simulate_run", _ok_sim)

        def _boom(*_a: object, **_k: object) -> object:
            msg = "same exception type as the simulate case"
            raise ValueError(msg)

        monkeypatch.setattr(stages, "measure_run", _boom)

        outcome = stages.run_point(self._ctx(tmp_path), backend="disk")
        assert outcome.status is RunStatus.MEASUREMENT_ERROR
        assert isinstance(outcome.exception, ValueError)

    def test_nonzero_exit_is_not_an_exception(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        import tidal.measurement._run_stages as stages

        def _no_probe(*_a: object, **_k: object) -> tuple[None, dict[str, object]]:
            return None, {}

        def _failing_sim(*_a: object, **_k: object) -> tuple[int, float, None]:
            return 3, 0.5, None

        monkeypatch.setattr(stages, "probe_for_run", _no_probe)
        monkeypatch.setattr(stages, "simulate_run", _failing_sim)

        outcome = stages.run_point(self._ctx(tmp_path), backend="disk")
        assert outcome.status is RunStatus.SIMULATION_FAILED
        assert outcome.exit_code == 3
        assert outcome.exception is None
