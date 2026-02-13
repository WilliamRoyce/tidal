# Hamiltonian from Lagrangian: Derivation and Implementation

Reference document for the energy measurement module's Hamiltonian
computation. This explains *why* the formula works, what assumptions it
makes, and where the subtleties lie — particularly around canonical vs.
simulation momenta.

## 1. The Problem

The TIDAL simulation pipeline evolves coupled PDEs derived from a
Lagrangian.  To measure energy conservation and wave conversion
probabilities, we need the Hamiltonian H — the conserved energy
functional.

**Naive approach (wrong for gauge theories):** Compute H = T + V by
manually writing kinetic + gradient + mass terms for each field.  This
fails for systems with derivative coupling (e.g., scalar-vector
coupling, Chern-Simons) because cross-field interaction energy lives in
gradient and cross-derivative operators, not in the identity-based
coupling matrix.

**Our approach:** Reconstruct H automatically from the Euler-Lagrange
equations stored in the JSON spec.  No manual per-system energy
formulas.

## 2. Canonical vs. Simulation Momenta

This is the critical conceptual subtlety.  **Read this carefully.**

### What the simulation stores

The simulation uses a first-order-in-time reformulation:

```
  d_t phi_i = pi_i            (definition)
  d_t pi_i  = RHS_i(phi, pi)  (from Euler-Lagrange equations)
```

Here `pi_i = d_t phi_i` for ALL fields — scalars, spatial vector
components, everything.  We call this the **simulation momentum**
`pi_sim`.

### What the Hamiltonian requires (in principle)

The Legendre transform H = sum_i pi_canonical_i * d_t phi_i - L
uses the **canonical momentum**:

```
  pi_canonical_i = dL / d(d_t phi_i)
```

For **scalar fields** (phi), this gives:

```
  pi_canonical = d_t phi = pi_sim    (identical!)
```

For **spatial vector components** (A_i) in a gauge theory, the
Lagrangian contains `F_{0i} = d_t A_i - d_i A_0`, and:

```
  pi_canonical_i = dL / d(d_t A_i) = F^{0i} = d_t A_i - d_i A_0
                 = pi_sim_i - d_i A_0    (NOT the same as pi_sim!)
```

The difference is `d_i A_0` — the spatial gradient of the temporal gauge
component.

### Why this doesn't matter for the Hamiltonian

This is the key insight.  When we compute H = sum pi_canonical * phi_dot - L
for the vector sector:

**Step 1:** Kinetic term from the Legendre transform.

```
  pi_canonical_i * d_t A_i = (d_t A_i - d_i A_0) * d_t A_i
                            = (d_t A_i)^2 - (d_i A_0)(d_t A_i)
```

**Step 2:** Kinetic term from -L.

```
  -L contains +1/2 * (d_t A_i - d_i A_0)^2
  = 1/2 * [(d_t A_i)^2 - 2(d_i A_0)(d_t A_i) + (d_i A_0)^2]
```

**Step 3:** Combine.

```
  H kinetic = pi_can * phi_dot - 1/2 F_{0i}^2
  = [(d_t A_i)^2 - (d_i A_0)(d_t A_i)]
    - 1/2[(d_t A_i)^2 - 2(d_i A_0)(d_t A_i) + (d_i A_0)^2]
  = 1/2 (d_t A_i)^2 - 1/2 (d_i A_0)^2
  = 1/2 pi_sim^2   - 1/2 |grad(A_0)|^2
```

The cross-terms `(d_i A_0)(d_t A_i)` cancel exactly!

### The result

The Hamiltonian uses **simulation momenta** (not canonical), plus a
correction for the constraint field's self-energy:

```
  H = 1/2 * sum_{dynamical} integral pi_sim^2 dV
    - 1/2 * sum_{constraint} integral |grad(C)|^2 dV
    + V_potential
```

The negative sign on the constraint gradient energy is the residue of
the canonical-vs-simulation momentum difference.  It arises physically
from the Minkowski metric signature g^{00} = -1: temporal components
have opposite-sign self-energy compared to spatial components.

**Bottom line:** `pi_sim = d_t phi` is the correct momentum for computing
the Hamiltonian kinetic energy.  The gauge correction (`d_i A_0`) is
NOT missing — it shows up as the negative constraint field self-energy.

## 3. The Virial Formula

### Euler's theorem for the potential

For quadratic Lagrangians (all systems in the TIDAL pipeline), the
potential energy V[phi] is a degree-2 functional of the fields.
Euler's homogeneous function theorem for functionals states:

```
  sum_i integral phi_i * (delta V / delta phi_i) dV = 2V
```

Since the Euler-Lagrange equations give:

```
  d^2_t phi_i = RHS_i = -(delta V / delta phi_i) + gyroscopic terms
```

we can write (excluding gyroscopic and velocity-dependent terms):

```
  V = -1/2 * sum_{i: dynamical} integral phi_i * RHS_i^{spatial} dV
```

This is what we call the **virial potential** `V_virial`.

### What the virial captures

The virial formula automatically includes:

1. **Self-gradient energy**: 1/2 |grad(phi)|^2
2. **Self-mass energy**: 1/2 m^2 phi^2
3. **Identity coupling**: g * phi * chi (cross-field mass coupling)
4. **Derivative coupling**: g * grad(phi) . A (scalar-vector coupling)
5. **Cross-derivative coupling**: any bilinear operator between fields
6. **Constraint-mediated coupling**: e.g., kCS * d_y(A_0) in A_1's eq

No manual enumeration of these terms is needed — the virial formula
reads them all from the RHS of the equations of motion.

### Proof for the single Klein-Gordon case

Consider `RHS = laplacian(phi) - m^2 * phi`.  Then:

```
  V_virial = -1/2 integral phi * [laplacian(phi) - m^2 * phi] dV
           = -1/2 integral phi * laplacian(phi) dV + 1/2 m^2 integral phi^2 dV
```

Integration by parts (or equivalently, Parseval's theorem with FFT):

```
  -integral phi * laplacian(phi) dV = integral |grad(phi)|^2 dV
```

So:

```
  V_virial = 1/2 integral |grad(phi)|^2 dV + 1/2 m^2 integral phi^2 dV
```

This is exactly the standard Klein-Gordon potential energy.

## 4. Constraint Field Self-Energy

### The sign flip

Constraint fields (time_derivative_order = 0) in the TIDAL pipeline
are always temporal gauge components (A_0 in electrodynamics, N or
N_i in gravity).  Their self-energy has the **opposite sign** compared
to scalar or spatial vector fields:

```
  V_constraint = sum_{j: constraint} [ -1/2 integral |grad(C_j)|^2 dV
                                       -1/2 m_j^2 integral C_j^2 dV ]
```

### Physical origin

The Minkowski metric g^{00} = -1 means that the "kinetic" term for A_0
in the Lagrangian has the wrong sign compared to spatial components:

```
  L contains: +1/2 g^{00} g^{ii} (d_i A_0)^2 = -1/2 (d_i A_0)^2
```

So the gradient energy of A_0 appears with a minus in L, and a minus
in H (since the A_0 equation has no d^2_t, it enters the Hamiltonian
purely as potential energy with the sign from -L).

For the mass term, the convention `mass_matrix[i][j] =
-(coefficient of identity(field_j) in equation_i)` means the A_0 row
stores `-Am2` (negative!) for a +Am2 mass-squared.  This correctly
gives a negative mass energy contribution for temporal components.

### Derivation from the Legendre transform

Starting from:

```
  L = 1/2 (d_t A_i - d_i A_0)^2 - 1/2 M^2 (A_0^2 - A_i^2) + ...
```

After the Legendre transform (Section 2), the A_0-dependent terms are:

```
  H contains: -1/2 |grad(A_0)|^2 - 1/2 M^2 A_0^2
```

Both are negative.  This is not an error — it reflects the indefinite
nature of the energy functional in gauge theories.  Total energy is
still positive because the constraint equation relates A_0 to the
spatial fields (e.g., via Gauss's law), ensuring that V_constraint +
V_virial > 0 for physical configurations.

## 5. Excluded Terms

### Gyroscopic forces (`first_derivative_t`)

Terms proportional to `d_t phi_j` in the equation for `phi_i` are
gyroscopic forces (velocity-dependent, configuration-independent).  They
do no work because they appear in antisymmetric pairs:

```
  A_1 equation: +kCS * d_t(A_2)   (Chern-Simons coupling)
  A_2 equation: -kCS * d_t(A_1)   (antisymmetric partner)
```

Their contribution to d_t H cancels:

```
  d_t H contains: pi_1 * kCS * pi_2 + pi_2 * (-kCS) * pi_1 = 0
```

For `first_derivative_t(A_0)` in the phi equation: A_0 is a
constraint field with zero momentum, so `d_t(A_0)` is determined by
the constraint equation's time evolution, not by a dynamical momentum.
This term contributes zero to the virial at each instant.

### Momentum references (`pi_N`)

Terms like `gradient_x(pi_0)` in the A_1 equation reference the
momentum of another field.  These are velocity-dependent forces, not
potential energy contributions.  For constraint fields, `pi_N = 0`
(zero momentum by construction).

## 6. Complete Formula

Putting it all together:

```
  H = 1/2 * sum_{i: dynamical} integral pi_sim_i^2 dV      [kinetic]
    + V_virial                                                [potential: virial]
    + V_constraint_self                                       [potential: constraint]
```

where:

```
  V_virial = -1/2 * sum_{i: dynamical} integral phi_i * RHS_i^{spatial} dV
```

with RHS_i^{spatial} excluding `first_derivative_t` and `pi_N` terms, and:

```
  V_constraint_self = sum_{j: constraint} [-1/2 integral |grad(C_j)|^2 dV
                                           -1/2 m_j^2 integral C_j^2 dV]
```

### Backward compatibility

For systems without constraints (pure scalars, coupled scalars):
- V_constraint_self = 0
- V_virial = sum of per-field gradient + mass + identity coupling
- Equivalent to the original formula

For systems with constraints (gauge theories):
- V_constraint_self adds the missing temporal component self-energy
- V_virial captures all derivative coupling automatically
- Energy conservation is restored

## 7. Assumptions and Limitations

1. **Quadratic Lagrangian**: The virial formula V = -1/2 sum phi * RHS
   is exact only for degree-2 potentials.  All current TIDAL systems
   are quadratic.  For higher-order Lagrangians, the virial theorem
   generalizes to V = -1/n * sum phi * RHS for degree-n, but this is
   not implemented.

2. **Constant coefficients**: Position-dependent mass or coupling terms
   (coordinate_dependent in JSON) cause ValueError.  The virial still
   holds mathematically, but the numerical integration would need to
   evaluate the position-dependent coefficient at each grid point.

3. **Minkowski-like signature**: The constraint sign flip assumes
   g^{00} = -1.  For Euclidean signature or non-standard conventions,
   the sign logic would need revision.

4. **No nonlinear terms**: The formula does not handle phi^3, phi^4, or
   other nonlinear potentials.  These would require explicit integration
   of the potential density, not the virial shortcut.

5. **Dirichlet BCs + cross_derivative: non-self-adjoint discretization**.
   The discrete `cross_derivative_xy` operator with Dirichlet ghost
   cells is NOT self-adjoint at boundary cells.  Specifically, the
   matrix element `M[(0,j),(0,j+1)] = +1/(4 dx dy)` while
   `M[(0,j+1),(0,j)] = -1/(4 dx dy)` — opposite signs.  This makes
   the discrete curl-curl system non-Hamiltonian, so no quadratic
   energy functional is exactly conserved.  The resulting energy drift
   (~30%) occurs even when using py-pde's own operators to compute the
   virial.  With periodic BCs the same system conserves energy to
   machine precision (~1e-10).  This is a fundamental discretization
   limitation, not a bug in the energy module.

   **Affected example:** `coupled_proca/` — two massive vector fields
   in a Dirichlet cavity.  The `measure_coupling.py` script uses a
   relaxed threshold (0.5) and documents this as an expected limitation.

   **Workarounds:**
   - Use periodic boundary conditions where physically appropriate
   - Use higher-resolution grids (concentrates the asymmetry to a
     thinner boundary layer)
   - Accept the energy drift as a discretization artifact and focus on
     conversion probability measurements, which remain valid

## 8. Verification Cases

### Case 1: Single Klein-Gordon field

```
  L = 1/2 (d_t phi)^2 - 1/2 |grad(phi)|^2 - 1/2 m^2 phi^2
  RHS = laplacian(phi) - m^2 phi
  V_virial = -1/2 phi * (laplacian phi - m^2 phi) = 1/2|grad phi|^2 + 1/2 m^2 phi^2
  V_constraint = 0
  H = 1/2 pi^2 + 1/2|grad phi|^2 + 1/2 m^2 phi^2    [standard result]
```

### Case 2: Coupled scalars (identity coupling)

```
  L = 1/2(d_t phi)^2 + 1/2(d_t chi)^2 - 1/2|grad phi|^2 - 1/2|grad chi|^2
      - 1/2 m^2(phi^2 + chi^2) - g phi chi
  V_virial = 1/2|grad phi|^2 + 1/2 m^2 phi^2 + g phi chi
           + 1/2|grad chi|^2 + 1/2 m^2 chi^2
  interaction = V_virial - (self phi) - (self chi) = g phi chi
  H = 1/2 pi_phi^2 + 1/2 pi_chi^2 + V_virial         [standard result]
```

### Case 3: Scalar-vector coupling (the motivating example)

```
  L = 1/2(d_t phi)^2 + 1/2(d_t A_i - d_i A_0)^2 - 1/4 F_{ij}^2
      - 1/2 m_phi^2 phi^2 - 1/2 M^2(A_0^2 - A_i^2)
      + gSV phi div(A) + kCS/2 eps A dA

  H = 1/2 pi_phi^2 + 1/2 pi_{A1}^2 + 1/2 pi_{A2}^2    [kinetic: pi_sim!]
    + V_virial     [gradient + mass + gSV coupling + kCS coupling + A_0 cross-terms]
    + (-1/2|grad A_0|^2 - 1/2 M^2 A_0^2)               [constraint self-energy]
```

## 9. Implementation Map

| Concept | Function | Location |
|---------|----------|----------|
| First derivative (FFT/FD) | `_first_derivative` | `_energy.py:89` |
| Second derivative (FFT/FD) | `_second_derivative` | `_energy.py:111` |
| Gradient energy density | `_gradient_energy_density` | `_energy.py:138` |
| Spatial operator dispatch | `_apply_spatial_operator` | `_energy.py:156` |
| Momentum field detection | `_is_momentum_field` | `_energy.py:211` |
| Term coefficient resolution | `_resolve_term_coefficient` | `_energy.py:216` |
| Term target resolution | `_resolve_term_target` | `_energy.py:237` |
| Single-field energy | `compute_field_energy` | `_energy.py:278` |
| Mass matrix resolution | `_resolve_mass_squared` | `_energy.py:332` |
| **Virial potential** | `_compute_virial_potential` | `_energy.py:369` |
| **Constraint self-energy** | `_compute_constraint_self_energy` | `_energy.py:422` |
| **System Hamiltonian** | `compute_system_energy` | `_energy.py:456` |
