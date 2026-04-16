# Supervisor Meeting — 17 April 2026

**Period**: 10 Apr (last meeting) to 17 Apr 2026
**Version**: v0.28.0 → v0.31.3 (82 commits, 21 issues closed)

---

## Summary

Three workstreams this week:

1. **Plasma Gertsenshtein** — completed `gertsenshtein_proca` with Raffelt-Stodolsky validation
2. **Dark photon** — understood the 276-run null result (eigenvalue structure); built plasma extension with 5 HPC sweep scripts
3. **Practice talk** — 13-section Beamer deck with speaker notes

---

## 1. Plasma Gertsenshtein (`gertsenshtein_proca`)

**Theory**: `examples/gertsenshtein_proca/theory.toml`

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}R - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} - \frac{m_A^2}{2}a_\mu a^\mu$$

Einstein-Maxwell with an effective photon mass $m_A^2$ representing the plasma frequency $\omega_p^2$. The mass is applied to the perturbation field $a_\mu$ (not the full background $A_\mu$), following Domcke et al. (2025).

### Raffelt-Stodolsky formula (1988)

$$P(g \to \gamma) = \sin^2(2\theta)\,\sin^2\!\Bigl(\frac{\Delta_{\mathrm{osc}}\,D}{2}\Bigr), \qquad \tan 2\theta = \frac{\kappa B_0}{|\Delta|}, \quad \Delta = -\frac{m_A^2}{2\omega}$$

- **Massless limit** ($m_A^2 = 0$): $P = \sin^2(\kappa B_0 D/2)$ — standard Gertsenshtein (1962)
- **Off-resonance** ($m_A^2 \gg \kappa B_0 \omega$): conversion suppressed as $P \propto (\kappa B_0/m_A^2)^2$
- **Resonance**: when $m_A^2 = m_T^2$ (photon mass matches dark photon mass), $A_{\mathrm{osc}}$ is maximised

The numerical solver reproduces this formula to < 0.04% across the full $m_A^2$ range.

### Literature

| Claim                            | Citation                                       | Reference                                                                                       |
| -------------------------------- | ---------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| Effective photon mass for plasma | Domcke, Garcia-Cely & Lee 2025                 | arXiv:2507.16609, Eq. (1484); motivation lines 115--117                                         |
| Mass on perturbation field       | Domcke et al. 2025                             | Eq. (187) defines $j_{\mathrm{eff}}^\mu$ on perturbation; Eq. (1484) adds $\mu^2$ on same field |
| Two-state Schrodinger mixing     | Berlin, Gonzalez-Solis, Melville, Trickle 2024 | arXiv:2405.08865, Eq. (309)                                                                     |
| Lorentzian resonance             | Raffelt & Stodolsky 1988                       | Phys. Rev. D 37, 1237                                                                           |

**Methodological note**: Domcke et al. add the mass at the EOM level. We add it at the Lagrangian level and derive the same EOM by variational calculus — equivalent for a quadratic mass term.

---

## 2. Dark Photon Kinetic Mixing — The Null Result

### Lagrangian (torsion trace as dark photon)

**Theory**: `examples/torsion_dark_photon/theory.toml`

$$\mathcal{L} = \frac{1}{\kappa^2}R + \alpha\, T^\lambda{}_{\lambda\mu}\, T^{\nu\mu}{}_\nu - \frac{\xi}{4}(F_T)_{\mu\nu}(F_T)^{\mu\nu} + \delta_m\, F_{\mu\nu}(F_T)^{\mu\nu} - \frac{1}{4}F_{\mu\nu}F^{\mu\nu}$$

The torsion trace vector $T_\mu \equiv T^\lambda{}_{\lambda\mu}$ acts as a dark photon with:

- Mass from $\alpha I_3$ — the proper torsion-trace mass invariant from the PGT action (not an ad hoc Proca term)
- Yang-Mills kinetic term $-\frac{\xi}{4}(F_T)^2$ where $(F_T)_{\mu\nu} = \partial_\mu T_\nu - \partial_\nu T_\mu$
- Kinetic mixing $\delta_m F \cdot F_T$ with the photon

Plain $R$ (not $\tilde{R}$): using $\tilde{R}$ (Riemann-Cartan) introduced an unavoidable torsion mass floor $m_T \sim 1/\kappa \gg \omega_G$ via the ChangeCurvature decomposition, pushing torsion deep into the adiabatic regime at all $B_0$. Switching to plain $R$ eliminates this.

### HPC campaign result

**Exact null**: $P_{\max}(h_5 \to a_1) = \sin^2(\kappa B_0 t/2)$ to $6.7 \times 10^{-6}$ across the entire $(\alpha, \xi, \delta_m)$ parameter space.

### Why the null is correct

This is **not** simply Holdom triviality. Holdom (1986) applies to massless fields; the dark photon here is massive. The real mechanism:

After diagonalising the kinetic matrix (which has off-diagonal entries from $\delta_m F \cdot F_T$), the mass eigenstates have asymmetric mixing. The off-diagonal coupling in the equations of motion — i.e. how strongly $a$ sources $t$ and vice versa through mass mixing — is:

$$M_{ta} = \frac{-2\,\delta_m\,m_A^2}{4\delta_m^2 - \xi}$$

where $m_A^2$ is the photon mass. In the vacuum ($m_A^2 = 0$), $M_{ta} = 0$ exactly. The photon is massless: after kinetic diagonalisation, the mass eigenstate for $a$ remains an exact zero-mass mode with no cross-coupling to $t$. The dark photon sits in an orthogonal eigenstate that the TT graviton initial condition cannot populate.

**This is an algebraic, not numerical, result** — it holds at all values of $\alpha$, $\xi$, $\delta_m$.

### What breaks this

A photon effective mass $m_A^2 \neq 0$ (plasma) makes $M_{ta} \neq 0$, rotating the eigenstates so the dark photon becomes accessible from the Gertsenshtein channel. This motivated the plasma extension.

---

## 3. Dark Photon Plasma Model

**Theory**: `examples/dark_photon_plasma/theory.toml`

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}R + \alpha_3\, T^\lambda{}_{\lambda\mu}\, T^{\nu\mu}{}_\nu - \frac{\xi}{4}(F_T)_{\mu\nu}(F_T)^{\mu\nu} + \delta_m\, F_{\mu\nu}(F_T)^{\mu\nu} - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} - \frac{m_A^2}{2}a_\mu a^\mu$$

Direct extension of the vacuum CDT model: adds $-\frac{m_A^2}{2}a_\mu a^\mu$ (photon plasma mass on the perturbation field) to break the eigenmode degeneracy.

### Parameters and sweep ranges

| Symbol     | Role                                         | Range        |
| ---------- | -------------------------------------------- | ------------ |
| $m_A^2$    | Photon plasma mass                           | 0.01 -- 5.0  |
| $\alpha_3$ | PGT torsion trace mass ($m_T^2 = 2\alpha_3$) | 0.005 -- 2.5 |
| $\delta_m$ | Kinetic mixing                               | -2.0 -- 2.0  |
| $\xi$      | Torsion trace kinetic coefficient            | 0.1 -- 5.0   |

Ghost-freedom: $|\delta_m| < \sqrt{\xi}/2$.

### Amplification metrics

$$A_{\mathrm{total}} = \frac{P_{\max}}{P_{\mathrm{GR}}}, \qquad A_{\mathrm{dark}} = \frac{P_{\max}(\delta_m \neq 0)}{P_{\max}(\delta_m = 0)}$$

$A_{\mathrm{dark}} > 1$ means the dark photon enhances total conversion above the plasma baseline; $< 1$ means it drains energy from the photon channel.

### References

- An, Pospelov, Pradler (2013), arXiv:1302.3884 — dark photon conversion in plasma
- Holdom (1986), Phys. Lett. B 166, 196 — kinetic mixing
- Raffelt & Stodolsky (1988), Phys. Rev. D 37, 1237 — two-state mixing framework

---

## Questions for Supervisor

### Dark photon model

- [ ] **Model extensions**: the vacuum null is clean — any thoughts on the most physically motivated next term to add? Our current thinking is this model stays as a pure dark photon sweep, separate from the PGT torsion geometry sweeps.

### PhD applications

- [ ] **Will Handley** PhD — when to email? Before or after the practice talk?
- [ ] **Sven Krippendorf**: no reply yet. You mentioned emailing on my behalf if no reply by tomorrow (17 Apr) — is that still the plan?

### Practice talk

- [ ] Any feedback on structure / emphasis before rehearsal?
