"""Gauge-orbit chains against the captured E.cal pencil (GH #477 stage F1).

Builds the linearized-diffeomorphism and U(1) orbit chains in the modal
solver's own basis (full-fftn, slot-major: index = slot*N + mode) and
measures the chain identities

    A·R_0 = 0,   A·R_{j+1} = B·R_j,   B·R_d = 0

whose defect is, by the linearized Noether identity, exactly linear in
the background-EOM violation (see the arc-F plan / GH #477).

Conventions verified against tidal/solver/modal.py:
  * full-fftn basis, ``y_hat[slot] = np.fft.fftn(field_profile).ravel()``
  * ``k = 2*pi*fftfreq(N, d=dx)``
  * E.cal spec coordinate "x" IS physical z (plane-wave reduction along z)

Component map (chart (t,x,y,z)): h_0=h_tt h_1=h_tx h_2=h_ty h_3=h_tz
h_4=h_xx h_5=h_xy h_6=h_xz h_7=h_yy h_8=h_yz h_9=h_zz; a_0..a_3=a_t..a_z.
Under the plane-wave reduction only d_t and d_z survive.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from tidal.solver import modal
from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_convolution_matrix_with_constraints,
    _build_k_axes_full,
    _build_k_grid,
    _pencil_deflate,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

REPO = Path(__file__).resolve().parents[3]
SPEC = REPO / "examples/data/gertsenshtein_ungauged_e_dual_gaussian.json"
# Frozen Phase-E geometry (scripts/hpc_submit_drafts/v3e_localised/_geometry.env)
BASE_PARAMS = {"kappa": 1.0, "Bpeak": 0.01, "sigB": 5.0, "zc1": 25.0, "zc2": 75.0}
LENGTH = 100.0


@dataclass
class Ctx:
    """Captured pencil plus everything needed to build orbit chains."""

    A: np.ndarray
    B: np.ndarray
    N: int
    z: np.ndarray
    k: np.ndarray
    slot: dict[str, int]
    n_slots: int
    params: dict[str, float]

    @property
    def dim(self) -> int:
        return self.n_slots * self.N

    def B0(self) -> np.ndarray:
        """Background magnetic field B_x(z) = -d_z Abar_y."""
        p = self.params
        g1 = np.exp(-((self.z - p["zc1"]) ** 2) / (2 * p["sigB"] ** 2))
        g2 = np.exp(-((self.z - p["zc2"]) ** 2) / (2 * p["sigB"] ** 2))
        return p["Bpeak"] * (g1 - g2)

    def Abar(self) -> np.ndarray:
        """Background potential Abar_y(z) (covariant component)."""
        from scipy.special import erf

        p = self.params
        pre = -p["Bpeak"] * p["sigB"] * np.sqrt(np.pi / 2)
        e1 = erf((self.z - p["zc1"]) / (np.sqrt(2) * p["sigB"]))
        e2 = erf((self.z - p["zc2"]) / (np.sqrt(2) * p["sigB"]))
        return pre * (e1 - e2)

    def pack(self, comps: dict[str, np.ndarray]) -> np.ndarray:
        """Real-space per-field profiles -> state vector in the fftn basis."""
        y = np.zeros(self.dim, dtype=np.complex128)
        for name, prof in comps.items():
            s = self.slot[name]
            y[s * self.N : (s + 1) * self.N] = np.fft.fftn(prof)
        return y


def capture(n: int, bpeak: float | None = None) -> Ctx:
    """Capture the raw (A, B) pencil the modal builder hands to deflation."""
    params = dict(BASE_PARAMS)
    if bpeak is not None:
        params["Bpeak"] = bpeak
    spec = EquationSystem.from_dict(json.loads(SPEC.read_text()))
    grid = GridInfo(bounds=((0.0, LENGTH),), shape=(n,), periodic=(True,))
    layout = StateLayout.from_spec(spec, n)
    ce = CoefficientEvaluator(spec, grid, params)
    k_grid = _build_k_grid(_build_k_axes_full(grid))

    cap: dict[str, np.ndarray] = {}

    def _capture(A, B, *, context="", diagnostics=None, tag=("pencil", 0)):  # noqa: ANN001
        cap["A"], cap["B"] = np.array(A, copy=True), np.array(B, copy=True)
        return np.asarray(A, dtype=np.complex128), np.eye(
            A.shape[0], dtype=np.complex128
        )

    modal._pencil_deflate = _capture  # noqa: SLF001
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, _, _, o2r, _ = _build_convolution_matrix_with_constraints(
                spec, layout, grid, ce, k_grid, tuple(grid.shape)
            )
    finally:
        modal._pencil_deflate = _pencil_deflate  # noqa: SLF001

    slot: dict[str, int] = {f: o2r[s] for f, s in layout.field_slot_map.items()}
    slot.update({"v_" + f: o2r[s] for f, s in layout.velocity_slot_map.items()})
    return Ctx(
        A=cap["A"],
        B=cap["B"],
        N=n,
        z=grid.axes_coords(0),
        k=np.asarray(2.0 * np.pi * np.fft.fftfreq(n, d=grid.dx[0])),
        slot=slot,
        n_slots=cap["A"].shape[0] // n,
        params=params,
    )


# --------------------------------------------------------------------------
# Orbit chains.  R_j holds the coefficient of the j-th time derivative of the
# gauge function: y(t) = sum_j R_j f^(j)(t).  A field's velocity slot picks up
# the term one level up (v_f = d_t f), which is the load-time chain rule the
# production consumer will use instead of exporting velocities.
# --------------------------------------------------------------------------

GENERATORS = ("xi_t", "xi_x", "xi_y", "xi_z", "chi_u1")


def _profile(ctx: Ctx, mode: int) -> tuple[np.ndarray, np.ndarray]:
    """Real-space gauge profile for Fourier mode ``mode`` and its d_z."""
    hat = np.zeros(ctx.N, dtype=np.complex128)
    hat[mode] = 1.0
    return np.fft.ifftn(hat), np.fft.ifftn(1j * ctx.k * hat)


def chain(
    ctx: Ctx, gen: str, mode: int, *, covariant: bool = False
) -> list[np.ndarray]:
    """Chain [R_0, ..., R_d] for one generator and one gauge-profile mode.

    ``covariant`` selects delta a_mu = xi^nu Fbar_{nu mu} (the representative
    with no compensating U(1) piece) over the full Lie derivative
    delta a_mu = L_xi Abar_mu.  The two differ by a U(1) transformation with
    chi = xi.Abar, and the Lorenz gauge-fixing term in this spec breaks them
    differently -- see GH #477 foundation 4.
    """
    xi, dxi = _profile(ctx, mode)
    levels: list[dict[str, np.ndarray]] = [{}, {}, {}]

    if gen == "xi_t":
        # dh_tt = 2 d_t xi_t ; dh_tz = d_z xi_t
        levels[0]["h_3"] = dxi
        levels[1]["h_0"] = 2.0 * xi
        levels[1]["v_h_3"] = dxi
    elif gen == "xi_x":
        # dh_tx = d_t xi_x ; dh_xz = d_z xi_x
        levels[0]["h_6"] = dxi
        levels[1]["h_1"] = xi
        levels[1]["v_h_6"] = dxi
    elif gen == "xi_y":
        # dh_ty = d_t xi_y ; dh_yz = d_z xi_y ; photon via Abar (y is the
        # only nonzero background component)
        levels[0]["h_8"] = dxi
        levels[1]["h_2"] = xi
        levels[1]["v_h_8"] = dxi
        if covariant:
            # xi^y Fbar_{y z} = -xi_y * Abar' = +xi_y * B0
            levels[0]["a_3"] = ctx.B0() * xi
            levels[1]["v_a_3"] = ctx.B0() * xi
        else:
            ab = ctx.Abar()
            levels[0]["a_3"] = ab * dxi
            levels[1]["v_a_3"] = ab * dxi
            levels[1]["a_0"] = ab * xi
            levels[2]["v_a_0"] = ab * xi
    elif gen == "xi_z":
        # dh_tz = d_t xi_z ; dh_zz = 2 d_z xi_z ; da_y = xi^z d_z Abar_y
        # (both representatives agree here: Abar_z = 0)
        levels[0]["h_9"] = 2.0 * dxi
        levels[0]["a_2"] = -ctx.B0() * xi
        levels[1]["h_3"] = xi
        levels[1]["v_h_9"] = 2.0 * dxi
        levels[1]["v_a_2"] = -ctx.B0() * xi
        levels[2]["v_h_3"] = xi
    elif gen == "chi_u1":
        # da_mu = d_mu chi  (a_1, a_2 untouched under the reduction)
        levels[0]["a_3"] = dxi
        levels[1]["a_0"] = xi
        levels[1]["v_a_3"] = dxi
        levels[2]["v_a_0"] = xi
    else:  # pragma: no cover - guarded by GENERATORS
        raise ValueError(f"unknown generator {gen!r}")

    while levels and not levels[-1]:
        levels.pop()
    return [ctx.pack(c) for c in levels]


def equilibrate(
    A: np.ndarray, B: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Row-equilibrate the pencil (the engine's own left preconditioner)."""
    s = np.maximum(np.max(np.abs(A), axis=1), np.max(np.abs(B), axis=1))
    s[s == 0.0] = 1.0
    return A / s[:, None], B / s[:, None], s


def defect(A: np.ndarray, B: np.ndarray, R: list[np.ndarray]) -> np.ndarray:
    """Stacked chain-identity defect: [A R_0; A R_{j+1} - B R_j; B R_d]."""
    parts = [A @ R[0]]
    parts.extend(A @ R[j + 1] - B @ R[j] for j in range(len(R) - 1))
    parts.append(B @ R[-1])
    return np.concatenate(parts)


def sigma(A: np.ndarray, B: np.ndarray, R: list[np.ndarray]) -> float:
    """Relative determination measure ||D(u)|| / ||C(u)|| (plan step 2d)."""
    c = np.linalg.norm(np.concatenate(R))
    return float(np.linalg.norm(defect(A, B, R)) / max(c, 1e-300))
