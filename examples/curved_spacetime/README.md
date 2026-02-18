# Curved Spacetime Examples

This directory contains examples of field equations on curved (non-flat) spacetimes.

## Phase 1: Static Conformal Factor

**Files**: `conformal_static.toml`

A Klein-Gordon scalar field on conformally flat spacetime with **constant** conformal factor.

### Physics

Metric: `ds^2 = Omega^2 * (-dt^2 + dx^2)` where `Omega = 2` (constant)

For constant `Omega`:
- Christoffel symbols = 0 (derivatives of constant metric vanish)
- Effective mass: `m_eff^2 = m^2 * Omega^2 = 1 * 4 = 4`
- Wave equation: `d2_t phi = d2_x phi - m_eff^2 phi`

### Verification

Static conformal with `Omega=2, m=1` should produce identical dynamics to flat Klein-Gordon with `m^2=4`.

### Run

```bash
# Derive equations
tidal derive examples/curved_spacetime/conformal_static.toml

# Run simulation (via run.sh which handles both conformal and de Sitter)
cd examples/curved_spacetime && bash run.sh
```

---

## Phase 2: de Sitter Expansion (Time-Dependent)

**Files**: `de_sitter.toml`

A Klein-Gordon scalar field on de Sitter spacetime with **time-dependent** conformal factor.

### Physics

Metric: `ds^2 = Omega(t)^2 * (-dt^2 + dx^2)` where `Omega(t) = e^{Ht}`

For time-dependent `Omega(t) = e^{Ht}`:
- Christoffel symbols are non-zero (derivatives of metric != 0)
- Hubble friction term: `-n*H * d_t phi` where n = number of spatial dimensions
- For 1+1D (n=1): `-H * d_t phi`
- For 2+1D (n=2): `-2H * d_t phi`
- Wave equation (1+1D): `d2_t phi = d2_x phi - H * d_t phi - m^2 * e^{2Ht} * phi`

The Hubble friction term causes wave damping during expansion:
- Amplitude decays as `exp(-H*t/2)` for 1+1D
- Energy decays as `exp(-H*t)` for 1+1D

### Key Features

- **Non-zero Christoffels**: Computed by xAct from the metric definition
- **First-order time derivative**: New `first_derivative_t` operator for Hubble friction
- **Dimension-dependent friction**: Coefficient is `-n*H` where n = spatial dimensions
- **Time-dependent coefficients**: Mass term scales with `e^{2Ht}`

### Run

```bash
# Derive equations
tidal derive examples/curved_spacetime/de_sitter.toml

# Run simulation (via run.sh which handles both conformal and de Sitter)
cd examples/curved_spacetime && bash run.sh
```
