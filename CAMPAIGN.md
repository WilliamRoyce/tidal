# PGT Survey Campaign

**Central question:** Can PGT torsion amplify the Gertsenshtein effect? If so, which
Lagrangian sector and at what parameters?

**Plan file:** `/home/vscode/.claude/plans/binary-snacking-chipmunk.md`
**Branch:** `hpc/pgt-survey` (created from `feature/torsion-sweeps`; rebase onto `main` once PR merges)

**Fixed params throughout:** `kappa=1.0`, `B0=0.01`
**IC convention:** `--ic plane-wave --ic-amplitude 1e-2 --periodic --bounds 0:50`
**Paired runs:** every theory gets an amplification run (`P_max:maximize`) + suppression run (`P_max:minimize`) submitted simultaneously.

---

## Theory Registry

| ID | Campaign name | JSON | Free params | IC fields |
|----|--------------|------|-------------|-----------|
| T1 | Dark-Photon-Plasma | `dark_photon_plasma.json` | mA2, deltam, xi, alpha3 | h=?, a=? (inspect) |
| T2 | Einstein-Cartan | `torsion_gertsenshtein_b5_zero.json` | alpha1, alpha2, alpha3 | h=?, a=? (inspect) |
| T3 | R²-PGT | `torsion_gertsenshtein.json` | alpha1-3, b5 | h=?, a=? (inspect) |
| T4 | Ricci-EM | `torsion_gertsenshtein_nonminimal.json` | alpha1-3, delta1 | h=?, a=? (inspect) |
| T5 | YM-PGT | `torsion_gertsenshtein_general_nonminimal.json` | beta1-3, xi, delta1, chi, zeta1-3 | h=?, a=? (inspect) |
| T6 | YM-PGT-CP | `torsion_gertsenshtein_parity_odd.json` | beta1-3, xi, delta1, chi, zeta1-3, d14-21, zt1-6 | h=?, a=? (inspect) |
| T7 | Complete-Even-PGT | `pgt_complete_even.json` (TBD, needs derivation) | T5 + chi2-10, b2-6 | TBD |
| T8 | Complete-Odd-PGT | `pgt_complete_odd.json` (TBD, needs derivation) | T6 + chitilde1-36, d1-13 | TBD |
| EH | Einstein-Maxwell-EH | TBD (GitHub issue filed) | c1, c2 | TBD |

**Sub-theories of T5 (no new JSON needed, fix unused params to 0):**
| Sub-ID | Campaign label | Free params | Fixed to 0 | Paper |
|--------|---------------|-------------|------------|-------|
| T5-Barker | Barker-PGT | beta1-3, xi, chi | delta1, zeta1-3 | Barker 2406.12826 |
| T5-Shapiro | Shapiro-PGT | beta1-3, zeta1-3 | xi, delta1, chi | Shapiro hep-th/0103093 |
| T5-Bahamonde | Bahamonde-PGT | beta1-3, delta1 | xi, chi, zeta1-3 | Bahamonde et al. |

---

## Stage Progress

### Stage 0: Local Smoke Tests (MANDATORY gate before any HPC submission)

| Theory | Inspected | IC fields found | Simulated | Energy OK | Notes |
|--------|-----------|----------------|-----------|-----------|-------|
| T1 Dark-Photon-Plasma | [ ] | h=?, a=? | [ ] | [ ] | |
| T2 Einstein-Cartan | [ ] | h=?, a=? | [ ] | [ ] | |
| T3 R²-PGT | [ ] | h=?, a=? | [ ] | [ ] | T3 extra: check b5 correction locally at 3 values |
| T4 Ricci-EM | [ ] | h=?, a=? | [ ] | [ ] | |
| T5 YM-PGT | [ ] | h=?, a=? | [ ] | [ ] | |
| T6 YM-PGT-CP | [ ] | h=?, a=? | [ ] | [ ] | |

### Stage A: Dark-Photon-Plasma 4D nested sampling

- [ ] Stage 0 gate passed for T1
- [ ] HPC amplification job submitted (job ID: ?)
- [ ] HPC suppression job submitted (job ID: ?)
- [ ] Results pulled
- [ ] Analysis: D_KL(params) computed — threshold > 0.05 nats
- Notes:

### Stage B: Einstein-Cartan null (T2)

- [ ] Stage 0 gate passed for T2
- [ ] HPC null job submitted (job ID: ?)
- [ ] Results pulled
- [ ] Analysis: D_KL < 0.005 nats (null confirmed) or signal found
- Notes:

### Stage C: R²-PGT perturbative b5 (T3)

- [ ] Stage 0: local b5 correction check (b5 ∈ {0, 1e-3, 5e-3}, P_max(b5≠0) ≠ P_max(b5=0))
- [ ] Stage C gate passed
- [ ] HPC order=0 amplification job (job ID: ?)
- [ ] HPC order=1 amplification job (job ID: ?)
- [ ] Results pulled
- [ ] Analysis: D_KL(b5, ord=1) − D_KL(b5, ord=0) > 0.01 nats
- Notes:

### Stage D: Nonminimal survey

#### D1: Ricci-EM / T4

- [ ] Stage 0 gate passed for T4
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- [ ] Analysis: D_KL(delta1) > 0.05 nats?
- Notes:

#### D2: YM-PGT / T5 + paper sub-theories

- [ ] Stage 0 gate passed for T5
- [ ] HPC: full YM-PGT amplify (job ID: ?)
- [ ] HPC: full YM-PGT suppress (job ID: ?)
- [ ] HPC: Barker-PGT amplify (job ID: ?)
- [ ] HPC: Shapiro-PGT amplify (job ID: ?)
- [ ] Results pulled
- [ ] Analysis: any D_KL > 0.1 nats?
- Notes:

#### D3: YM-PGT-CP / T6

- [ ] Stage 0 gate passed for T6
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- [ ] Analysis: any D_KL > 0.1 nats?
- Notes:

### Stage E: New theory derivations

#### E1: Complete-Even-PGT / T7

- [ ] theory.toml designed (`examples/pgt_complete_even/`)
- [ ] Wolfram derivation complete
- [ ] Smoke test: non-zero chi2+ or b2+ correction confirmed
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- Trigger: Stage D shows signal in χ or R×R direction

#### E2: Complete-Odd-PGT / T8

- [ ] theory.toml designed (`examples/pgt_complete_odd/`)
- [ ] Wolfram derivation complete
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- Trigger: Stage D3 shows parity-odd signal

---

## HPC Job Log

| Job ID | Theory | Run type | Status | Wall time | Notes |
|--------|--------|----------|--------|-----------|-------|
| — | — | — | — | — | — |

---

## Key Findings

*(Filled in as results arrive)*

- Stage A: 
- Stage B: 
- Stage C: 
- Stage D1: 
- Stage D2: 
- Stage D3: 
- Stage E: 

---

## GitHub Issues Filed

| Issue | Topic |
|-------|-------|
| #296 | Euler-Heisenberg F⁴ extension (#271 resolved, defer EH theory to future) |
| #297 | Non-propagating/constraint torsion theories (single irreducible sector investigation) |
