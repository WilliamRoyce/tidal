# v3_permissive campaign submit scripts

**Architecture:** v3 (post-2026-05-08 supervisor pivot) — see [`docs/V3_ARCHITECTURE.md`](../../../docs/V3_ARCHITECTURE.md).

**Common settings across all scripts**:

- **Likelihood:** `P_max:maximize` (amp) / `P_max:minimize` (sup). Sim runs every sample. No probe gate (--gated NOT set), no Hwang-Noh, no upper P_max cap.
- **Soft floor:** `logL = -100 + Normal(0, 1.0)` for sim divergence / NaN / exception. Tunable via `--soft-floor-noise SIGMA`.
- **Compactified priors:** sign-symmetric dimensionless couplings → `arctan_uniform:-89:89`; positive-definite kinetic → `log_uniform:1e-3:1e3`.
- **Geometry:** plane-wave + uniform B₀ as in v2 (Phase E will replace with localised geometry).
- **Probe metadata:** γ_eff, k_dominant, n_tachyonic_modes, borderline_stability all flow through chain CSV for post-hoc analysis (no longer gates).

## Script inventory

| Script | Theory | Mode | Wall budget |
|---|---|---|---|
| `d1_amp.sh` | D1 Ricci-EM | `P_max:maximize` | 8 h std |
| `d1_sup.sh` | D1 Ricci-EM | `P_max:minimize` | 12 h std |
| `stage_a_amp.sh` | Dark-Photon-Plasma | `P_max:maximize` | INTR feasible |
| `stage_a_sup.sh` | Dark-Photon-Plasma | `P_max:minimize` | INTR feasible |
| `d20_bahamonde_amp.sh` | T5 Bahamonde sub | `P_max:maximize` | INTR / 2 h std |
| `d20_bahamonde_sup.sh` | T5 Bahamonde sub | `P_max:minimize` | INTR / 2 h std |
| `d21_barker_amp.sh` | T5 Barker sub | `P_max:maximize` | INTR / 2 h std |
| `d21_barker_sup.sh` | T5 Barker sub | `P_max:minimize` | INTR / 2 h std |
| `d22_shapiro_amp.sh` | T5 Shapiro sub | `P_max:maximize` | INTR / 3 h std |
| `d22_shapiro_sup.sh` | T5 Shapiro sub | `P_max:minimize` | INTR / 3 h std |
| `d23_full_amp.sh` | T5 full 9-D | `P_max:maximize` | INTR (1 session) |
| `d23_full_sup.sh` | T5 full 9-D | `P_max:minimize` | INTR (1 session) |

All D2 chains use INTR (starts immediately vs days of standard queue wait).

## Policy: 1 INTR session per chain (updated 2026-05-15)

**Each chain = exactly 1 INTR session (≤1 h). Accept whatever PolyChord state exists at TIMEOUT.**
This is a landscape-overview survey, not a convergence-required result.

- If a chain accumulates < 5K dead points by TIMEOUT (gross under-sampling), allow **1 resume**.
- Otherwise: pull → plot → compare → submit-next. No further resumes.
- **No `*_resume.sh` scripts going forward.**

Rationale: v3 sup chains are extremely multi-modal (401 clusters for 5-param D2.0 sup after 7 sessions)
and the PolyChord precision criterion is impractical to reach in finite INTR budget. The landscape
overview — where in parameter space suppression is possible — is the scientific deliverable.

Observed v3 timings (updated 2026-05-15):

- D2.0 amp (5p): 22 min INTR ✓ (18K dead pts, ESS=8338) — 1 session
- D2.0 sup (5p): **1-session truncated** at 58.5K dead pts, 174/401 active clusters (7-session outlier; new policy prevents repetition)
- D2.1 amp (6p): 27 min INTR ✓ (21K dead pts, ESS=9517, 27 clusters) — 1 session
- D2.1 sup (6p): ~60 min → submit 1 session, take result
- D2.2 amp/sup (8p): ~60 min → submit 1 session each
- D2.3 amp/sup (9p): ~60 min → submit 1 session each

If/when the cubed-sphere joint prior scripts land, sibling `scripts/hpc_submit_drafts/v3_jointprior/` will contain `--joint-prior` versions for direct comparison.

## Usage

Each script is `set -euo pipefail` and self-contained. Submit via:

```bash
bash scripts/hpc_submit_drafts/v3_permissive/d1_amp.sh
```

Pre-flight: ensure HPC venv is up to date (`bash scripts/hpc_refresh_venv_tar.sh` after any `uv sync` or refactor that changes `_likelihood.py`).
