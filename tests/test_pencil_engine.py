"""Analytic validation of the ordered-QZ pencil engine (GH #457 Stage 3).

``_pencil_deflate`` is the general DAE engine for the modal solver: it
must reproduce hand-derived reduced dynamics on synthetic pencils BEFORE
being wired into any builder. Cases cover the structure classes the
corpus exhibits and the ones it does not yet (generality is the point —
the engine has no per-case formulas to get wrong):

* trivial and invertible-B pencils (must equal B⁻¹A exactly);
* an index-2 constraint chain shaped like the E-class (a promoted row
  with no own-velocity content whose velocity is determined through a
  mixed row) — fully hand-solved reduced generator;
* a #260-style rank-deficient mass matrix (position elimination through
  K) — hand-solved effective frequency;
* a genuinely singular (gauge) pencil — must refuse loudly, never pick
  an evolution silently.
"""

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver.modal import (
    SingularPencilError,
    _pencil_deflate,  # pyright: ignore[reportPrivateUsage]
)


def _rand_complex(rng: np.random.Generator, n: int) -> np.ndarray:
    return rng.standard_normal(n) + 1j * rng.standard_normal(n)


class TestTrivialPencils:
    def test_identity_b_returns_a(self) -> None:
        rng = np.random.default_rng(3)
        A = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
        A_eff, proj = _pencil_deflate(A, np.eye(4, dtype=np.complex128))
        np.testing.assert_allclose(A_eff, A, atol=1e-12)
        np.testing.assert_allclose(proj, np.eye(4), atol=1e-12)

    def test_invertible_diagonal_b(self) -> None:
        rng = np.random.default_rng(4)
        A = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        B = np.diag(np.array([2.0, -0.5, 1.5], dtype=np.complex128))
        A_eff, proj = _pencil_deflate(A, B)
        np.testing.assert_allclose(A_eff, np.linalg.solve(B, A), atol=1e-12)
        np.testing.assert_allclose(proj, np.eye(3), atol=1e-12)


class TestIndexTwoChain:
    """E-class-shaped analytic pencil, fully hand-solved.

    Slots y = (qφ, vφ, qc, vc)::

        q̇φ = vφ
        v̇φ = −ω²·qφ + γ·vc          (mixed row references the promoted
                                      field's velocity)
        q̇c = vc
        0   = a·qc + b·vφ + d·qφ     (promoted row: algebraic LHS, no own
                                      velocity — D̃_cc = 0, index 2)

    Hand solution: qc = −(d·qφ + b·vφ)/a and, closing the algebraic loop
    through v̇φ,  vc = (b·ω²·qφ − d·vφ)/(a + γ·b), giving the reduced
    on-manifold oscillator::

        v̇φ = −ω²·a/(a+γb) · qφ  −  γ·d/(a+γb) · vφ
    """

    OMEGA2 = 1.7
    A_C = 0.9
    B_C = 0.35
    D_C = -0.6
    GAMMA = 0.4j  # complex, like an i·k spatial coupling

    def _pencil(self) -> tuple[np.ndarray, np.ndarray]:
        w2, a, b, d, g = self.OMEGA2, self.A_C, self.B_C, self.D_C, self.GAMMA
        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        B[0, 0] = 1.0
        A[0, 1] = 1.0
        B[1, 1] = 1.0
        A[1, 0] = -w2
        A[1, 3] = g
        B[2, 2] = 1.0
        A[2, 3] = 1.0
        A[3, 2] = a
        A[3, 1] = b
        A[3, 0] = d
        return A, B

    def _manifold_state(self, qphi: complex, vphi: complex) -> np.ndarray:
        w2, a, b, d, g = self.OMEGA2, self.A_C, self.B_C, self.D_C, self.GAMMA
        qc = -(d * qphi + b * vphi) / a
        vc = (b * w2 * qphi - d * vphi) / (a + g * b)
        return np.array([qphi, vphi, qc, vc], dtype=np.complex128)

    def test_manifold_dimension_is_two(self) -> None:
        A, B = self._pencil()
        _A_eff, proj = _pencil_deflate(A, B)
        assert round(float(np.real(np.trace(proj)))) == 2

    def test_reduced_dynamics_match_hand_solution(self) -> None:
        w2, a, b, d, g = self.OMEGA2, self.A_C, self.B_C, self.D_C, self.GAMMA
        A, B = self._pencil()
        A_eff, proj = _pencil_deflate(A, B)
        rng = np.random.default_rng(457)
        for _ in range(4):
            qphi, vphi = _rand_complex(rng, 2)
            y = self._manifold_state(qphi, vphi)
            # On-manifold states are fixed by the projector.
            np.testing.assert_allclose(proj @ y, y, atol=1e-10)
            dy = A_eff @ y
            vdot_phi = (-w2 * a / (a + g * b)) * qphi + (-g * d / (a + g * b)) * vphi
            np.testing.assert_allclose(dy[0], vphi, atol=1e-10)
            np.testing.assert_allclose(dy[1], vdot_phi, atol=1e-10)
            np.testing.assert_allclose(dy[2], y[3], atol=1e-10)
            # v̇c follows by differentiating the vc recovery.
            vdot_c = (b * w2 * vphi - d * vdot_phi) / (a + g * b)
            np.testing.assert_allclose(dy[3], vdot_c, atol=1e-10)

    def test_flow_stays_on_manifold(self) -> None:
        from scipy.linalg import expm

        A, B = self._pencil()
        A_eff, _proj = _pencil_deflate(A, B)
        y0 = self._manifold_state(0.8, -0.3)
        y_t = expm(A_eff * 2.5) @ y0
        # The evolved state still satisfies the algebraic row exactly.
        a, b, d = self.A_C, self.B_C, self.D_C
        residual = a * y_t[2] + b * y_t[1] + d * y_t[0]
        assert abs(residual) < 1e-9


class TestRankDeficientMass:
    """#260-style: M = [[1,1],[1,1]] (rank 1), K = diag(−k1², −k2²).

    Hand solution: constraint −k1²q1 + k2²q2 = 0 (row difference), so
    q2 = (k1²/k2²)·q1, and the surviving oscillator has
    ω² = k1²k2²/(k1² + k2²).
    """

    K1SQ = 2.0
    K2SQ = 0.5

    def _pencil(self) -> tuple[np.ndarray, np.ndarray]:
        # slot order: q1, v1, q2, v2
        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        B[0, 0] = 1.0
        A[0, 1] = 1.0
        B[2, 2] = 1.0
        A[2, 3] = 1.0
        # M rows on the velocity slots
        B[1, 1] = 1.0
        B[1, 3] = 1.0
        A[1, 0] = -self.K1SQ
        B[3, 1] = 1.0
        B[3, 3] = 1.0
        A[3, 2] = -self.K2SQ
        return A, B

    def test_effective_frequency(self) -> None:
        A, B = self._pencil()
        A_eff, proj = _pencil_deflate(A, B)
        c = self.K1SQ / self.K2SQ
        omega2 = self.K1SQ * self.K2SQ / (self.K1SQ + self.K2SQ)
        assert round(float(np.real(np.trace(proj)))) == 2
        rng = np.random.default_rng(260)
        for _ in range(3):
            (q1, v1) = _rand_complex(rng, 2)
            y = np.array([q1, v1, c * q1, c * v1], dtype=np.complex128)
            np.testing.assert_allclose(proj @ y, y, atol=1e-10)
            dy = A_eff @ y
            np.testing.assert_allclose(dy[0], v1, atol=1e-10)
            np.testing.assert_allclose(dy[1], -omega2 * q1, atol=1e-10)
            np.testing.assert_allclose(dy[2], c * v1, atol=1e-10)
            np.testing.assert_allclose(dy[3], -c * omega2 * q1, atol=1e-10)


class TestSingularPencilGaugeCompletion:
    """Singular pencils get an EXPLICIT, warned gauge completion.

    Dependent equations are dropped and the least-constrained state
    combinations frozen at IC — the loud form of the min-norm/frozen
    convention the pre-#457 machinery applied silently (needed because
    sweeps legitimately cross #260-class critical points, and GH #465
    specs ship redundant equation copies). SingularPencilError is
    reserved for pencils the completion cannot regularize.
    """

    def test_disconnected_pair_frozen_with_warning(self) -> None:
        # Slots (qφ, vφ, qc, vc): a clean oscillator plus a (qc, vc)
        # pair whose equation row is missing entirely — vc undetermined.
        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        B[0, 0] = 1.0
        A[0, 1] = 1.0
        B[1, 1] = 1.0
        A[1, 0] = -1.0
        B[2, 2] = 1.0
        A[2, 3] = 1.0  # kinematic row of the orphan pair
        with pytest.warns(UserWarning, match="quotiented out"):
            A_eff, _proj = _pencil_deflate(A, B, context="unit test")
        # The oscillator dynamics are untouched by the completion.
        y = np.array([0.7, -0.2, 0.0, 0.0], dtype=np.complex128)
        dy = A_eff @ y
        np.testing.assert_allclose(dy[0], -0.2, atol=1e-10)
        np.testing.assert_allclose(dy[1], -0.7, atol=1e-10)
        assert SingularPencilError is not None  # class stays exported

    def test_pinned_directions_are_recorded_per_slot(self) -> None:
        """GH #468 pin + certify: the quotient reports WHAT it pinned.

        The orphan (qc, vc) pair is the Kronecker chain x(λ) = (0,0,1,λ);
        both of its slots are pinned, the oscillator's slots are not — so
        a measurement on (qφ, vφ) is certifiable and one on (qc, vc) is
        flagged. The per-slot overlap is the squared row norm of the
        recorded orthonormal pinned basis.
        """
        from tidal.solver.modal import PencilDiagnostics

        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        B[0, 0] = 1.0
        A[0, 1] = 1.0
        B[1, 1] = 1.0
        A[1, 0] = -1.0
        B[2, 2] = 1.0
        A[2, 3] = 1.0
        diag = PencilDiagnostics()
        with pytest.warns(UserWarning, match="quotiented out"):
            _pencil_deflate(
                A, B, context="unit test", diagnostics=diag, tag=("evolution", 0)
            )
        assert len(diag.entries) == 1
        (tag, pinned, tau) = diag.entries[0]
        assert tag == ("evolution", 0)
        assert tau > 0
        assert pinned.shape == (4, 2)
        row_overlap = np.sum(np.abs(pinned) ** 2, axis=1)
        np.testing.assert_allclose(row_overlap[:2], 0.0, atol=1e-12)
        np.testing.assert_allclose(row_overlap[2:], 1.0, atol=1e-12)
