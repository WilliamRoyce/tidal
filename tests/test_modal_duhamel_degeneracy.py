"""Near-degeneracy sweep for the v6 Pass 1 Duhamel kernel.

Validates ``tidal.solver.modal._duhamel_kernel`` across the full
single-to-far-from-degenerate range of ``|(μ − λ)·t|``. The closed-form
Duhamel integral loses accuracy in the direct ``expm1(z)/z`` formula as
|z| → 0 (catastrophic cancellation) and in the Taylor fallback as
|z| → 1 (truncation at order 12). The module's threshold of 1e-5 is
chosen so both branches remain within machine precision of their
analytic answer.

References
----------
Al-Mohy & Higham 2011, SIAM J. Sci. Comput. 33, §3 Table 3.1.
"""

# ruff: noqa: RUF002, ERA001 — math symbols in docstrings; commented stubs.

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver.modal import _duhamel_kernel


def _reference_kernel(lam: complex, mu: complex, t: float) -> complex:
    """High-precision reference using mpmath's mpf/expm1 (50 digits).

    The reference evaluates the exact closed-form
    ``G(λ, μ; t) = [exp(μt) − exp(λt)] / (μ − λ)`` at 50-digit
    precision. For μ = λ it returns the degenerate limit ``t·exp(λt)``
    without going through the singular expression.
    """
    try:
        from mpmath import (
            exp as mp_exp,  # type: ignore[import-untyped]
        )
        from mpmath import mp, mpc  # type: ignore[import-untyped]
    except ImportError:  # pragma: no cover — defer to numpy fallback
        if mu == lam:
            return complex(t * np.exp(lam * t))
        return complex((np.exp(mu * t) - np.exp(lam * t)) / (mu - lam))

    mp.dps = 50
    lam_m = mpc(lam.real, lam.imag)
    mu_m = mpc(mu.real, mu.imag)
    t_m = mp.mpf(t)

    if mu_m == lam_m:
        val = t_m * mp_exp(lam_m * t_m)
    else:
        val = (mp_exp(mu_m * t_m) - mp_exp(lam_m * t_m)) / (mu_m - lam_m)
    return complex(float(val.real), float(val.imag))


class TestDuhamelKernelScalar:
    """Scalar-input edge cases."""

    def test_exact_degeneracy(self) -> None:
        """G(λ, λ; t) = t·exp(λt) exactly (no 0/0)."""
        lam = 0.5 + 1.3j
        t = 2.0
        expected = t * np.exp(lam * t)
        got = _duhamel_kernel(lam, lam, t)
        assert abs(got - expected) < 1e-14, f"error {abs(got - expected):.2e}"

    def test_zero_t(self) -> None:
        """G(λ, μ; 0) = 0 for any λ, μ."""
        for lam, mu in [(0.0, 0.0), (1j, 2j), (-1.0 + 0.3j, 0.5 - 0.7j)]:
            got = _duhamel_kernel(lam, mu, 0.0)
            assert abs(got) < 1e-14, f"G(·, ·, 0) should be 0, got {got}"

    def test_purely_imaginary_far_apart(self) -> None:
        """Harmonic oscillator eigenvalues: match analytic formula."""
        lam = 1j
        mu = 3j
        t = 1.0
        expected = (np.exp(mu * t) - np.exp(lam * t)) / (mu - lam)
        got = _duhamel_kernel(lam, mu, t)
        assert abs(got - expected) < 1e-14

    def test_vectorised_broadcast(self) -> None:
        """Array lam/mu broadcast elementwise."""
        lam = np.array([1j, 2j, 3j], dtype=np.complex128)
        mu = lam.copy()  # exact degeneracy on diagonal
        lam_i = lam[:, None]
        mu_j = mu[None, :]
        t = 1.5
        G = _duhamel_kernel(lam_i, mu_j, t)
        # Diagonal: t·exp(λt)
        for i in range(3):
            expected_diag = t * np.exp(lam[i] * t)
            assert abs(G[i, i] - expected_diag) < 1e-14, (
                f"diagonal[{i}]: {G[i, i]} vs {expected_diag}"
            )
        # Off-diagonal: direct formula
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                expected = (np.exp(mu[j] * t) - np.exp(lam[i] * t)) / (mu[j] - lam[i])
                assert abs(G[i, j] - expected) < 1e-13, (
                    f"[{i},{j}]: {G[i, j]} vs {expected}"
                )


class TestDuhamelKernelDegeneracySweep:
    """|μ − λ|·t swept from 1e-15 to 1.0; error ≤ 1e-13 throughout."""

    @pytest.mark.parametrize(
        "log10_dz",
        [-15.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.0, -1.0, 0.0],
    )
    def test_real_axis(self, log10_dz: float) -> None:
        """Real eigenvalues across the degeneracy region."""
        lam = 0.3  # arbitrary real baseline
        dz = 10**log10_dz
        mu = lam + dz
        t = 1.0
        got = _duhamel_kernel(lam, mu, t)
        expected = _reference_kernel(complex(lam, 0), complex(mu, 0), t)
        rel_err = abs(got - expected) / max(abs(expected), 1e-30)
        assert rel_err < 1e-13, (
            f"log10|z|={log10_dz}: got={got}, expected={expected}, "
            f"rel_err={rel_err:.3e}"
        )

    @pytest.mark.parametrize(
        "log10_dz",
        [-15.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.0, -1.0, 0.0],
    )
    def test_imaginary_axis(self, log10_dz: float) -> None:
        """Purely imaginary eigenvalues (Hamiltonian spectra)."""
        lam = 1j  # baseline oscillator
        dz = 10**log10_dz
        mu = lam + dz * 1j
        t = 1.0
        got = _duhamel_kernel(lam, mu, t)
        expected = _reference_kernel(lam, mu, t)
        rel_err = abs(got - expected) / max(abs(expected), 1e-30)
        assert rel_err < 1e-13, (
            f"log10|z|={log10_dz}: got={got}, expected={expected}, "
            f"rel_err={rel_err:.3e}"
        )

    @pytest.mark.parametrize(
        "log10_dz",
        [-15.0, -10.0, -6.0, -3.0, 0.0],
    )
    @pytest.mark.parametrize("t", [0.1, 1.0, 10.0])
    def test_t_sweep(self, log10_dz: float, t: float) -> None:
        """Varying t preserves accuracy at fixed |μ−λ|."""
        lam = 0.5 + 1.3j
        dz = 10**log10_dz
        mu = lam + dz * (1.0 + 0.5j)
        got = _duhamel_kernel(lam, mu, t)
        expected = _reference_kernel(lam, mu, t)
        rel_err = abs(got - expected) / max(abs(expected), 1e-30)
        assert rel_err < 1e-13, (
            f"log10|dz|={log10_dz}, t={t}: got={got}, expected={expected}, "
            f"rel_err={rel_err:.3e}"
        )

    def test_crossover_threshold_continuity(self) -> None:
        """Direct and Taylor branches agree to < 1e-12 at the crossover."""
        from tidal.solver.modal import _PHI1_DEGENERACY_THRESHOLD

        lam = 1.0 + 0.5j
        t = 1.0
        threshold = _PHI1_DEGENERACY_THRESHOLD

        # Sample just below and just above the threshold; both should be
        # within 1e-12 of the 50-digit reference.
        for factor in (0.5, 0.9, 1.1, 2.0):
            dz = threshold * factor
            mu = lam + dz
            got = _duhamel_kernel(lam, mu, t)
            expected = _reference_kernel(lam, mu, t)
            rel_err = abs(got - expected) / abs(expected)
            assert rel_err < 1e-12, (
                f"factor={factor} (|z|≈{dz}): got={got}, "
                f"expected={expected}, rel_err={rel_err:.3e}"
            )
