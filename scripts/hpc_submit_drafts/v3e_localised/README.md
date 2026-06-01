# Phase E — Localized wavepacket + dual-Gaussian B-field submit scripts

Companion to `docs/PHASE_E_TRACKER.md` (primary tracker) and
`docs/V3_PHASE_E_DESIGN.md` (design rationale).

## Layout

- `_geometry.env` — **single source of truth** for geometry parameters. Every script
  here `source`s this file. Changing it triggers a re-baseline of all prior Phase E
  runs. See its header for the parameter rationale.
- `<stage>_<theory>_{amp,sup}.sh` — corner-plot chain submits (one pair per theory).
- `atlas/<stage>_<theory>/face_*.sh` — cubed-sphere atlas tile submits.
- `run_interactive_batch_*.sh` — parallel-launch wrappers for INTR interactive sessions.

## Workflow

1. Run `python scripts/v3e_boccaletti_preflight.py` — confirms `_geometry.env` is in
   the perturbative, away-from-node regime. Must exit 0 before any submit.
2. Push to HPC: `bash scripts/hpc_shuttle.sh push`.
3. Submit a single theory's corner pair: `bash run_interactive_batch_corner.sh e_cal`.
4. Or atlas: `bash run_interactive_batch_atlas.sh e_t2_einstein_cartan`.
5. After completion: `bash scripts/hpc_shuttle.sh pull <jobid>` writes the per-job
   README skeleton; analyst fills in verdict lines.

## Key choices (changing requires explicit re-baseline)

| Param | Value | Why |
|-------|-------|-----|
| `BPEAK` | 0.01 | κ·Bpeak·σ_B·√(2π)/2 ≈ 0.063 — perturbative, far from Boccaletti node |
| `R_ATLAS` | 0.4 | Brackets ~35% of T4 v3 MAP magnitude; perturbative regime |
| `K_CARRIER` | 3 | k·σ_B ≈ 15 → high-frequency Gertsenshtein |
| `T_END` | 60 | Wavepacket front clears B-field at t≈50; margin for trailing edge |
| `BC` | periodic | Required for modal solver (post-#367 fix in v0.42.0/0.42.1) |
| `SNAPSHOTS` | 2 | Mandatory per `feedback_snapshots_mandatory.md` |
