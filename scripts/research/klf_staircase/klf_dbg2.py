"""Sub-step equivalence bisection for preduce2/preduce4 on random pencils.

Entry precondition (post-preduceBF window): N = [[0_{n x m}, E_{n x n}],
[0_{p x m}, 0_{p x n}]] with E upper-triangular invertible; M arbitrary.
Every kernel must preserve Q' A0 Z = M, Q' B0 Z = N exactly at exit.
"""

from __future__ import annotations

import klf_port as kp
import numpy as np

rng = np.random.default_rng(7)


def make_window(n: int, m: int, p: int, complex_data: bool = True):
    rows, cols = n + p, m + n
    assert rows == cols, "square test"

    def rnd(r, c):
        x = rng.standard_normal((r, c))
        if complex_data:
            x = x + 1j * rng.standard_normal((r, c))  # noqa: PLR6104 — += breaks dtype promotion
        return x

    M = rnd(rows, cols).astype(kp.C128)
    N = np.zeros((rows, cols), dtype=kp.C128)
    E = np.triu(rnd(n, n)) + 3.0 * np.eye(n)
    N[:n, m:] = E
    return M, N


def check(label, A0, B0, M, N, Q, Z):
    em = np.linalg.norm(Q.conj().T @ A0 @ Z - M)
    en = np.linalg.norm(Q.conj().T @ B0 @ Z - N)
    print(f"  {label}: equiv M {em:.2e}  N {en:.2e}")
    return em + en


def run(kernel: str, n: int, m: int, p: int):
    M, N = make_window(n, m, p)
    A0, B0 = M.copy(), N.copy()
    Q = np.eye(M.shape[0], dtype=kp.C128)
    Z = np.eye(M.shape[1], dtype=kp.C128)
    tol = 1e-10
    if kernel == "p4":
        rho = kp.preduce4(n, m, p, M, N, Q, Z, tol, 0, 0, 0, 0)
        print(f"preduce4 n={n} m={m} p={p}: rho={rho}")
    elif kernel == "p2":
        tau, rho = kp.preduce2(n, m, p, M, N, Q, Z, tol, 0, 0, 0, 0)
        print(f"preduce2 n={n} m={m} p={p}: tau={tau} rho={rho}")
    elif kernel == "p3":
        # precondition for preduce3: same standard form
        rho = kp.preduce3(n, m, M, N, Q, Z, tol, 0, 0, 0, 0)
        print(f"preduce3 n={n} m={m}: rho={rho}")
    check("exit", A0, B0, M, N, Q, Z)
    # structural exit checks for p2: D block zero on first m-tau cols
    return M, N, Q, Z


if __name__ == "__main__":
    print("--- preduce4 alone (tau=0 path) ---")
    run("p4", 6, 3, 3)
    run("p4", 8, 2, 2)
    print("--- preduce3 alone ---")
    M, N = make_window(6, 3, 3)
    A0, B0 = M.copy(), N.copy()
    Q = np.eye(9, dtype=kp.C128)
    Z = np.eye(9, dtype=kp.C128)
    rho = kp.preduce3(6, 3, M, N, Q, Z, 1e-10, 0, 0, 3, 0)
    print(f"preduce3 n=6 m=3 (rtrail=3): rho={rho}")
    check("exit", A0, B0, M, N, Q, Z)
    print("--- preduce2 full (m>0) ---")
    run("p2", 6, 3, 3)
    run("p2", 10, 4, 4)
