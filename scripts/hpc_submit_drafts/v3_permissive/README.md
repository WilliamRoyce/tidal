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
| `d23_full_amp.sh` | T5 full 9-D | `P_max:maximize` | 6 h std |
| `d23_full_sup.sh` | T5 full 9-D | `P_max:minimize` | 6 h std |

Scripts default to `polychord_intr.sbatch` for the mid-res variants (smoke / cross-check) and `polychord_standard.sbatch` for hi-res publication runs. The publication-run scripts are committed; INTR variants can be derived locally by changing `--template` and reducing `--nlive`.

If/when the cubed-sphere joint prior (parallel session) lands, sibling `scripts/hpc_submit_drafts/v3_jointprior/` will contain `--joint-prior` versions of the same campaigns for direct comparison.

## Usage

Each script is `set -euo pipefail` and self-contained. Submit via:

```bash
bash scripts/hpc_submit_drafts/v3_permissive/d1_amp.sh
```

Pre-flight: ensure HPC venv is up to date (`bash scripts/hpc_refresh_venv_tar.sh` after any `uv sync` or refactor that changes `_likelihood.py`).
