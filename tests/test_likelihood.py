"""Inline Hwang-Noh perturbativity gate (``P_MAX_LIMIT``) tests.

The gate rejects samples whose simulation produces ``P_max > 0.5``
under a ``P_max:maximize`` or ``P_max:minimize`` likelihood — those
values fall outside the linearised regime (Hwang & Noh 2002, PRD 65,
124010).  Gaussian / threshold likelihoods may legitimately probe
parameter regions where ``P_max`` exceeds the limit, so the gate is
applied only to the maximise/minimise modes.

The gate condition is on the *raw* ``P_max``, even when the campaign
likelihood maximises the amplification ratio
``A = P_max / P_GR`` via ``--baseline-formula``: the Hwang-Noh limit
is on the absolute conversion probability, not on its ratio to
baseline.  See ``docs/tex/stability_probe.tex`` for the role of this
gate in the chain-time rejection cascade.
"""

from __future__ import annotations

from tidal.inference._likelihood import P_MAX_LIMIT, _is_non_perturbative

# Sample 391's P_max value from the Stage C audit (the slow-growth
# Hwang–Noh case the probe alone misses).  Used here as a concrete
# above-limit value, not as a regression on the audit itself.
SAMPLE_391_P_MAX = 9.82


class TestPerturbativityGate:
    def test_p_max_limit_constant(self) -> None:
        """Sanity: the perturbativity limit is the documented Hwang-Noh value."""
        assert P_MAX_LIMIT == 0.5

    def test_maximize_above_limit_rejected(self) -> None:
        """``P_max:maximize`` with sample 391's P_max ≫ 0.5 → reject."""
        assert _is_non_perturbative("P_max", "maximize", SAMPLE_391_P_MAX)

    def test_minimize_above_limit_rejected(self) -> None:
        """``P_max:minimize`` with P_max above limit → reject."""
        assert _is_non_perturbative("P_max", "minimize", 0.6)

    def test_maximize_at_limit_not_rejected(self) -> None:
        """At the boundary the linearisation is still considered valid."""
        assert not _is_non_perturbative("P_max", "maximize", P_MAX_LIMIT)

    def test_maximize_below_limit_not_rejected(self) -> None:
        assert not _is_non_perturbative("P_max", "maximize", 0.1)

    def test_gaussian_not_gated(self) -> None:
        """Gaussian likelihoods may fit P_max > 0.5; gate must not fire."""
        assert not _is_non_perturbative("P_max", "gaussian", SAMPLE_391_P_MAX)

    def test_threshold_not_gated(self) -> None:
        """Threshold likelihoods may exceed the perturbativity limit."""
        assert not _is_non_perturbative("P_max", "threshold", 1.0)

    def test_other_metric_not_gated(self) -> None:
        """The gate only applies to ``P_max``."""
        assert not _is_non_perturbative("L_mix", "maximize", SAMPLE_391_P_MAX)
