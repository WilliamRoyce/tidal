# Supervisor Meeting — 16 April 2026

**Period**: 10 Apr (last meeting) to 16 Apr 2026
**Version**: v0.28.0 → v0.31.3 (82 commits, 21 issues closed)

---

## Summary

Three workstreams this week:

1. **Plasma Gertsenshtein** — completed the `gertsenshtein_proca` example with Raffelt-Stodolsky validation
2. **Dark photon model** — explained the 276-run null result (eigenvalue structure, not just Holdom); built the plasma extension (`dark_photon_plasma`) with 5 HPC sweep scripts
3. **Modal solver unification** (#256--#260) — fixed rank-deficient kinetic matrices that were blocking the dark photon campaign

Also: practice talk slides written (13-section Beamer deck with speaker notes).

---

## 1. Plasma Gertsenshtein (completed example)

**Theory**: `examples/gertsenshtein_proca/theory.toml`

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}R - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} - \frac{m_A^2}{2}a_\mu a^\mu$$

Einstein-Maxwell with an effective photon mass $m_A^2$ representing plasma dispersion ($\omega_p^2 = m_A^2$). The mass term is on the **perturbation** field $a_\mu$ (not the full field $A_\mu$), matching the convention in Domcke et al. (2025).

### Analytical formula (Raffelt-Stodolsky 1988)

$$P(g \to \gamma) = \sin^2(2\theta)\,\sin^2\!\Bigl(\frac{\Delta_{\mathrm{osc}}\,D}{2}\Bigr)$$

where $\tan(2\theta) = \kappa B_0 / |\Delta|$, $\;\Delta = -m_A^2/(2\omega)$, $\;\Delta_{\mathrm{osc}} = \sqrt{\Delta^2 + (\kappa B_0)^2}$.

- **Massless limit** ($m_A^2 = 0$): reduces to $P = \sin^2(\kappa B_0 D/2)$ (Gertsenshtein 1962)
- **Large mass** ($m_A^2 \gg \kappa B_0 \omega$): conversion suppressed as $P \propto (\kappa B_0 / m_A^2)^2$

### Literature attribution

| Claim | Primary citation | Exact reference |
|-------|-----------------|-----------------|
| Effective photon mass models plasma | Domcke, Garcia-Cely & Lee 2025 | arXiv:2507.16609, Eq. (1484), "Towards Including Medium Effects" (line 1477); motivation lines 115--117 |
| Mass on perturbation field, not background | Domcke et al. 2025 | Eq. (187) defines $j_{\mathrm{eff}}^\mu$ on perturbation; Eq. (1484) adds $\mu^2$ on same field |
| Two-state Schrodinger mixing structure | Berlin, Gonzalez-Solis, Melville, Trickle 2024 | arXiv:2405.08865, Eq. (309), "Schrodinger-like equation" at line 311 |
| Lorentzian resonance validation | Raffelt & Stodolsky 1988 | Phys. Rev. D 37, 1237; cited in Domcke Eq. (612) line 618 and Berlin Eq. (335) lines 332--338 |

**Key methodological point**: Domcke et al. add the mass at the EOM level ($[\Box - \mu^2]A_h^\mu = -j_{\mathrm{eff}}^\mu$). We add it at the Lagrangian level ($-m_A^2/2\;a_\mu a^\mu$ with $a$ the xPert perturbation field) and derive the same EOM via variational calculus. Formally equivalent for a quadratic mass term.

**Status**: Theory derived, JSON validates, sweep scripts (B_0 family, mA2 detuning, resonance) tested and committed.

---

## 2. Dark Photon Kinetic Mixing — The Null Result

### What we tested (supervisor's suggestion, 2 Apr meeting)

Torsion trace vector as a dark photon with kinetic mixing $\delta_m F \cdot F_t$.

### 276-run HPC campaign result

**Exact null**: $P_{\max}(h_5 \to a_1) = \sin^2(\kappa B_0 t/2)$ to $6.7 \times 10^{-6}$ relative precision across the entire $(\alpha, \xi, \delta_m)$ parameter space.

### Why the null is correct

This is **not** simply Holdom triviality (which applies to massless fields). The dark photon IS massive (Proca mass $m_T^2$), so Holdom's field-redefinition argument doesn't directly apply. The real mechanism is:

1. **Asymmetric eigenvalue structure**: After diagonalising the kinetic matrix, the mass eigenstates have asymmetric mixing angles. The graviton-photon Gertsenshtein channel ($h_5 \leftrightarrow a_1$) projects onto mass eigenstates that are *orthogonal* to the dark photon eigenstate.

2. **Algebraic inaccessibility**: The $h_5$ initial condition (TT graviton) excites only the Gertsenshtein subspace. The dark photon subspace is algebraically decoupled --- not perturbatively suppressed, but exactly zero coupling from TT initial conditions.

3. **Eigenvalue picture**: The 3-field system ($h_5, a_1, t_1$) has eigenstates where the $h_5 \leftrightarrow a_1$ oscillation is an exact eigenmode. The dark photon $t_1$ sits in an orthogonal eigenstate that the Gertsenshtein IC cannot populate.

This is an exact statement about eigenvalue structure, not a perturbative or numerical one. Verified across 276 runs spanning 5 orders of magnitude in each parameter.

### What breaks this

A **photon effective mass** $m_A^2 \neq 0$ (from plasma) lifts the eigenmode degeneracy. When the photon and dark photon have *different* masses, the eigenstates rotate and kinetic mixing becomes genuinely observable:

$$E[t,a] = \frac{-2\,\delta_m\,m_A^2}{4\delta_m^2 - \xi} \neq 0 \quad\text{when } m_A^2 > 0$$

This motivated the plasma extension (Section 3).

---

## 3. Dark Photon Plasma Model (new this week)

**Theory**: `examples/dark_photon_plasma/theory.toml`

### Lagrangian

$$\mathcal{L} = \frac{1}{\kappa^2}R - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} - \frac{m_A^2}{2}a_\mu a^\mu - \frac{\xi}{4}F^{(t)}_{\mu\nu}F^{(t)\mu\nu} - \frac{m_T^2}{2}t_\mu t^\mu + \delta_m\,F_{\mu\nu}F^{(t)\mu\nu}$$

This is a **pure BSM dark photon** model --- no geometric/PGT interpretation. A fundamental vector $T_\mu$ with kinetic mixing to EM, in a plasma background ($m_A^2$ = photon effective mass).

### Parameters

| Symbol | Role | Sweep range |
|--------|------|-------------|
| $m_A^2$ | Photon effective mass (plasma proxy) | 0.01 -- 5.0 |
| $m_T^2$ | Dark photon Proca mass | 0.01 -- 5.0 |
| $\delta_m$ | Kinetic mixing strength | -2.0 -- 2.0 |
| $\xi$ | Dark photon kinetic coefficient | 0.1 -- 5.0 |
| $B_0$ | Background magnetic field | 0.001 -- 0.015 |

Ghost-freedom condition: $|\delta_m| < \sqrt{\xi}/2$.

### Amplification metrics

- $A_{\mathrm{total}} = P_{\max} / P_{\mathrm{GR}}$ where $P_{\mathrm{GR}} = \sin^2(\kappa B_0 t/2)$
- $A_{\mathrm{dark}} = P_{\mathrm{full}}(\delta_m \neq 0) / P_{\mathrm{plasma}}(\delta_m = 0)$ — isolates the dark photon contribution

$A_{\mathrm{dark}} > 1$ means dark photon enhances conversion; $< 1$ means it drains energy from the photon channel.

### HPC sweep scripts (5 configurations)

| Script | Type | Grid | Purpose |
|--------|------|------|---------|
| `sweep_B0.sh` | 1D, 20 points | $B_0 \in [0.001, 0.015]$ | B-field independence check |
| `sweep_mc.sh` | MC, 1000 LHS samples | 4D $(m_A^2, m_T^2, \delta_m, \xi)$ | Full parameter space exploration |
| `sweep_paired.sh` | MC baseline | Same grid, $\delta_m = 0$ | Plasma-only baseline for $A_{\mathrm{dark}}$ |
| `sweep_heatmap_mass.sh` | 2D, $50 \times 50$ | $m_A^2 \times m_T^2$ | Resonance structure visualisation |
| `sweep_heatmap_dxi.sh` | 2D, $50 \times 50$ | $\delta_m \times \xi$ | Coupling-kinetics landscape |

All scripts tested locally, ready for CSD3 submission.

### Key physics questions for the sweep

- Is there any $(m_A^2, m_T^2, \delta_m, \xi)$ region where $A_{\mathrm{dark}} > 1$?
- Does the Raffelt-Stodolsky resonance ($m_A^2 \approx m_T^2$) produce enhancement or just redistribution?
- Is the dominant effect suppression (energy drained into dark photon channel)?

### References

- An, Pospelov, Pradler (2013), arXiv:1302.3884 --- dark photon conversion formulas
- Holdom (1986), Phys. Lett. B 166, 196 --- kinetic mixing triviality
- Pospelov (2008), arXiv:0811.1030 --- dark photon portal Lagrangian
- Fabbrichesi et al. (2020), arXiv:2005.01515 --- dark photon review

---

## 4. Modal Solver Unification (#256--#260)

Rank-deficient kinetic matrices from torsion/Proca theories caused the modal solver to inject spurious tachyonic modes. This was producing factor-20 discrepancies between the CDT and FV dark photon formulations and blocking the campaign.

**Fixes applied** (5 issues, v0.31.0--v0.31.3):

- **#256**: Unified builder replacing 3 separate code paths (constraint, standard, generalized)
- **#257**: SVD-based Schur elimination for rank-deficient blocks
- **#258**: Generalized eigenvalue for asymmetric $M$ (reads `kinetic_coefficient_symbolic`, not hardcoded 1.0)
- **#259**: Robust eigenvector inverse via condition-number check
- **#260**: Pseudoinverse for rank-deficient $K_{cc}$ in constraint Schur complement

All 1,721 tests passing. Dark photon CDT and FV formulations now agree at all generic parameter points.

---

## 5. Practice Talk

- 13-section Beamer deck with UC Rev CMYK theme
- Speaker notes with 6 physics anecdotes (tachyon masquerade, constrained Hamiltonian, polarisation block-diagonal, volume element subtlety, Holdom null, propagating-torsion dead end)
- Files: `docs/talks/practice_2026_04/`
- **Status**: slides and notes written, needs rehearsal + supervisor feedback on structure/emphasis

---

## Questions for Supervisor

### Dark photon model

- [ ] The vacuum kinetic mixing null is clean and physically well-understood. Worth writing up as a negative result, or fold into a larger paper?
- [ ] For the plasma dark photon model: any thoughts on additional terms to include? We're treating this as a pure BSM dark photon (separate from PGT geometric sweeps).
- [ ] Should we run the plasma sweeps on CSD3 before or after the practice talk?

### PhD applications

- [ ] Will Handley PhD --- when should I email? Before or after the practice talk?
- [ ] Sven Krippendorf: no reply yet. You said you'd email on my behalf if no reply by tomorrow (17 Apr) --- is that still the plan?

### Practice talk

- [ ] Any feedback on slide structure / emphasis before I rehearse?

---

## Appendix: Commits by Theme

### Solver fixes (16 commits)
`18dfa9a` perf: sequential per-pair TraceBasisDummy (#246) |
`a66f30c` fix: epsilon index notation in parity-odd TOML (#245) |
`bf4395f` fix: defer torsion-based derived fields |
`7a9cb49` fix: NumericQ guard for symbolic selfCoeff |
`b708be9` fix: exclude algebraicDepFields from zero rules |
`9dac0e4` fix: revert depField exclusion |
`51a36fd` fix: bypass ToCanonical for deferred torsion (#255) |
`9f496ef` fix: restore h_5 kinetic terms via full-expand (#250) |
`7da64c7` fix: restrict LagrangianFullExpand to non-torsion |
`4fef016` refactor: unify deferred-field canonicalization |
`14e2b09` fix(modal): pre-solve B_lhs in generalized builder (#256) |
`273190e` refactor(modal): delete constraint path (#256) |
`1a8d870` fix(modal): generalized eig for rank-deficient M (#256) |
`8ddc9cb` fix(modal): rank-deficient asymmetric M (#258, #259) |
`7ed05e9` fix(modal): constraint recovery in SVD Schur (#260) |
`27351a8` fix(modal): pseudoinverse for rank-deficient K_cc (#260)

### Dark photon & plasma models (14 commits)
`8582161` feat: R-tilde to plain R in dark photon, add delta1 |
`cfaa655` refactor: TorsionCDT trace to fundamental vector |
`e48fd4d`--`6fe133e` Revert fundamental-vector (back to CDT after #256 fix) |
`93c56b4` feat: fundamental-vector dark photon as reference |
`1708f4b` feat: plasma dark-photon model with observable conversion |
`b032ac9` refactor: tighten MC sweep bounds + paired baseline |
`5a747dd` feat: proper sweep scripts with xi + amplification metrics |
Various sweep script updates and HPC recipes

### Documentation & talk (12 commits)
`8e6a6de` docs: pre-arXiv references |
`e7f6d7d` docs: CHANGELOG.md |
`858c3c7` docs: unified canonicalization and plasma physics |
`7080f57` examples: reproduce_figures.sh for gertsenshtein_proca |
`04ad3bd` docs: comprehensive dark photon torsion documentation (#262) |
`79517f0` docs: UC Rev CMYK theme for practice talk |
`5896101` docs: initial slides and references |
`b46c675` docs: modal solver update for #256 |
`69854f5` docs: rank-deficient pathology in troubleshooting |
`b806313` docs: CHANGELOG v0.31.0 |
`a824dea` fix: re-derive propagating model with per-pair TBD |
`a28a393` docs: derivation timing after per-pair TBD fix

### Infrastructure (10 commits)
`dead4bb` feat: unified plotting with publication-quality upgrades |
`43c6d81` fix: plotting integration + 18 new tests |
`79c97a5` feat: sweep-grouped plot type |
`4671627` feat: CSD3 shuttle script + sbatch templates |
`658d231` perf: pre-import heavy modules in worker init (7--18x HPC speedup) |
`d249031` chore: un-ignore theory JSONs |
Various HPC shuttle fixes

### Version bumps (13 commits)
v0.28.0 → v0.28.1 → v0.29.0 → v0.29.1 → v0.29.2 → v0.29.3 → v0.29.4 → v0.30.0 → v0.30.1 → v0.30.2 → v0.31.0 → v0.31.1 → v0.31.2 → v0.31.3

### Issues closed (21)
#122, #125, #126, #164, #176, #226, #229, #249, #250, #251, #252, #253, #254, #255, #256, #257, #258, #259, #260, #261, #262

(14 created and closed this week: #249--#262; 7 older issues also closed)
