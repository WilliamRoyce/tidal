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
| `d23_full_amp.sh` | T5 full 9-D | `P_max:maximize` | INTR / 1-2 resumes |
| `d23_full_sup.sh` | T5 full 9-D | `P_max:minimize` | INTR / 4-5 resumes |

All D2 chains use INTR (starts immediately vs days of standard queue wait). Amp chains
typically fit in one 1h INTR session; sup chains (finding the null region in wide v3 priors)
take 5-6× longer than v2 and will need successive resumes via `--read-resume`.

## INTR + resume workflow

If a chain times out (SLURM state TIMEOUT), resume it immediately:

```bash
# Hardcode the previous job's output directory
PREV_OUTPUT=/rds/user/wr286/hpc-work/tidal/hpc_results/<JOBID>/<chain_name>

bash scripts/hpc_shuttle.sh submit \
  --template scripts/hpc_templates/polychord_intr.sbatch \
  --name <same_name_as_before> --ntasks 76 --time 01:00:00 \
  --cmd "tidal sample examples/data/<json> \
    <same flags as original> \
    --output ${PREV_OUTPUT} \
    --read-resume"
```

PolyChord writes `_chains/tidal.resume` on every checkpoint; the new job picks up from there.
Repeat until sacct shows COMPLETED (not TIMEOUT). Each successive job writes to the SAME
output directory — pull once at the end with the original jobid as the src path.

Observed v3 timings (updated 2026-05-12):

- D2.0 amp (5p): 22 min INTR ✓ (18K dead pts, ESS=8338)
- D2.0 sup (5p): **5-8+ INTR sessions** — 190+ clusters after 2 sessions (43.7K dead pts), σ_logZ still growing; multi-modal discovery accelerating
- D2.1 amp (6p): 27 min INTR ✓ (21K dead pts, ESS=9517, 27 clusters)
- D2.1 sup (6p): ~90-180 min → 1-2 resumes (6 params + wide priors)
- D2.2 amp (8p): v2 was 38 min; v3 uncertain — plan for 1-2 resumes
- D2.2 sup (8p): ~150-300 min → 2-4 resumes
- D2.3 amp (9p): v2 was 1:09; v3 uncertain — plan for 1-3 resumes
- D2.3 sup (9p): **≥5 resumes** — 9-param sup + wide v3 priors = most compute-intensive chain

**Key finding**: v3 sup chains are extremely multi-modal (190+ PolyChord clusters for 5-param D2.0 theory).
The ×5 v2 slowdown estimate is now a **severe lower bound** — sup chains can take 5-8× longer.
Plan for 5-8+ INTR sessions per sup chain, 1-3 per amp chain.

If/when the cubed-sphere joint prior (parallel session) lands, sibling `scripts/hpc_submit_drafts/v3_jointprior/` will contain `--joint-prior` versions of the same campaigns for direct comparison.

## Usage

Each script is `set -euo pipefail` and self-contained. Submit via:

```bash
bash scripts/hpc_submit_drafts/v3_permissive/d1_amp.sh
```

Pre-flight: ensure HPC venv is up to date (`bash scripts/hpc_refresh_venv_tar.sh` after any `uv sync` or refactor that changes `_likelihood.py`).
