"""Versioned stability-probe profiles.

A ``StabilityProfile`` couples the three scientific decisions that
together determine whether a parameter point is rejected by the
pre-flight tachyonic guard:

* ``ic_mode``   — how the probe IC is constructed.
* ``threshold`` — the ``γ_eff`` cutoff in s⁻¹ above which a mode is
  classified tachyonic.
* ``t_test``    — the probe time over which ``γ_eff`` is averaged.

These three are coupled (e.g. tightening ``threshold`` without
matching ``ic_mode`` over-rejects high-k modes; see #323), so they
travel together as a named, immutable record. The profile name is
written to per-sample metadata so any chain CSV can be unambiguously
identified with the probe used to generate it.

Why a registry rather than free-form kwargs
-------------------------------------------

The Stage A v5 / D1 / D2 chains were all evaluated under
``(unit, 0.3, t_end)`` and cross-stage ``D_KL`` comparability requires
identifying that combination as one decision, not three. Allowing
callers to pass arbitrary kwarg combinations re-introduces the silent
regression that produced commit 25cc531 → 0581fdb (a threshold
tightening that broke 100 % of N=64 priors because the IC mode was
unchanged).

References
----------
- Issue #323 — threshold history, false-negative regime documentation.
- Commit 0581fdb — revert of 25cc531; root-cause analysis.
- ``tidal/measurement/_stability.py`` — probe implementation.
"""

from __future__ import annotations

import dataclasses
from typing import Literal


@dataclasses.dataclass(frozen=True)
class StabilityProfile:
    """A named, immutable probe configuration.

    Attributes
    ----------
    name : str
        Stable identifier persisted into chain metadata. Format
        ``v{N}-{ic_mode}-{threshold}`` by convention; canonical entries
        live in :data:`PROFILES`.
    ic_mode : {"unit", "consistent"}
        How the probe IC is constructed:

        * ``"unit"``       — ``y₀[src] = 1`` at the source slot,
          remaining slots zero. Legacy behaviour, matches commits prior
          to 0581fdb.
        * ``"consistent"`` — full plane-wave IC at ``ic_wavevector``
          passed through ``ensure_consistent_ic`` so constraint-field
          slots are populated exactly as the simulation sees them.
          High-k modes get ≈ 0 IC weight, eliminating the high-k
          false positives that forced commit 25cc531's revert.
    threshold : float
        ``γ_eff`` cutoff in s⁻¹. Mode is tachyonic if
        ``log(‖y(t_test)‖/‖y₀‖) / t_test > threshold``.
    t_test : float
        Probe time in seconds for the Padé evolution. ``γ_eff``'s
        dependence on ``t_test`` is second-order (it converges to
        ``max Re(λ)`` in the linear-eigenvalue limit), so callers
        usually pass the simulation's ``t_end`` to keep the probe
        consistent with what the run will see.
    perturbativity_p_max : float
        Hwang–Noh perturbativity bound used by the post-hoc audit's
        hard-reject gate. Theory-specific (graviton-photon: 0.5;
        likely tighter for parity-odd or higher-symmetry channels).
        Probe itself does not use this — it lives on the profile so
        audit policy travels with probe policy.
    notes : str
        Free-form provenance documentation. Surfaced in the registry
        listing and in error messages.
    """

    name: str
    ic_mode: Literal["unit", "consistent"]
    threshold: float
    t_test: float
    perturbativity_p_max: float = 0.5
    notes: str = ""

    @property
    def borderline_low(self) -> float:
        """Lower bound of the borderline strip used by audit selection.

        Empirically the false-negative regime is the top ~30 % of the
        ``[0, threshold]`` interval (see #323): ``γ_eff`` just below the
        threshold passes the probe but evolves to ``P_max > 0.5`` in
        the simulation. Audit prioritises this strip for re-simulation.
        """
        return 0.7 * self.threshold


PROFILES: dict[str, StabilityProfile] = {
    "v1-unit-0.3": StabilityProfile(
        name="v1-unit-0.3",
        ic_mode="unit",
        threshold=0.3,
        t_test=10.0,
        perturbativity_p_max=0.5,
        notes=(
            "Legacy probe; pinned for v5 + D1 + D2 + D3 cross-stage "
            "D_KL comparability. False-negative regime "
            "γ_eff ∈ (0.21, 0.30) at IC k documented in #323; "
            "mitigated by post-hoc audit (Stage B). "
            "DO NOT change while a campaign is in flight."
        ),
    ),
    "v2-consistent-0.1": StabilityProfile(
        name="v2-consistent-0.1",
        ic_mode="consistent",
        threshold=0.1,
        t_test=10.0,
        perturbativity_p_max=0.5,
        notes=(
            "Constraint-solved IC (matches simulation's "
            "ensure_consistent_ic); high-k modes get ≈ 0 IC weight so "
            "threshold can drop to 0.1 without N=64 over-rejection. "
            "NOT the campaign default until validation Stage C "
            "(simulation truth-table on disagreement set) completes "
            "and supervisor signs off. See #323 resolution path."
        ),
    ),
}


DEFAULT_PROFILE_NAME = "v1-unit-0.3"


def get_profile(name: str) -> StabilityProfile:
    """Look up a profile by name; raise ``KeyError`` with available list."""
    if name not in PROFILES:
        msg = (
            f"Unknown stability profile {name!r}. "
            f"Available: {sorted(PROFILES)}. "
            f"To add a new profile, edit "
            f"tidal/measurement/_stability_profile.py."
        )
        raise KeyError(msg)
    return PROFILES[name]


def adhoc_profile(
    *,
    threshold: float,
    t_test: float = 10.0,
    ic_mode: Literal["unit", "consistent"] = "unit",
    perturbativity_p_max: float = 0.5,
) -> StabilityProfile:
    """Construct an ad-hoc profile for tests / one-offs.

    Use only for unit tests and bisecting investigations. Production
    runs MUST reference a registered :data:`PROFILES` entry so the
    chain CSV's ``stability_profile`` column is meaningful across
    sessions.
    """
    name = f"adhoc-{ic_mode}-{threshold:g}-{t_test:g}"
    return StabilityProfile(
        name=name,
        ic_mode=ic_mode,
        threshold=threshold,
        t_test=t_test,
        perturbativity_p_max=perturbativity_p_max,
        notes="Ad-hoc profile (test or bisect). Do not commit chain results.",
    )
