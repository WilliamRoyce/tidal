# Dark Photon Torsion — HPC Sweep Recipes

Uses `scripts/hpc_shuttle.sh` and the generic
`scripts/hpc_templates/sapphire_cpu.sbatch` template. No dark-photon-
specific sbatch file is needed — the template handles it via the
`TIDAL_PARALLEL=112` env override baked into each shell command below.

All three sweep scripts are self-contained: they call `tidal sweep`,
then `tidal plot`, then `analyze_sweep.py`, and write all output to
`examples/data/torsion_dark_photon_sweep_<name>/`.

## Prerequisites

```bash
# One-time (or when local code changes)
scripts/hpc_shuttle.sh push      # rsync repo to CSD3 /rds/.../tidal
scripts/hpc_shuttle.sh setup     # create/update remote venv
```

## Sweep 1 — ξ scan (smallest; good smoke test)

20 points in ξ, ~20 runs × ~1s each. Runs comfortably on a single
sapphire node with `--qos=intr` (queue-free, 1-hour cap).

```bash
scripts/hpc_shuttle.sh submit \
    --template scripts/hpc_templates/intr_test.sbatch \
    --name dp_xi \
    --ntasks 112 \
    --time 00:30:00 \
    --cmd 'TIDAL_PARALLEL=112 bash examples/torsion_dark_photon/sweep_xi.sh'
```

Acceptance: 20/20 stable, `P_max ≈ 0.061208` for every point (ξ does
not affect the h_5→a_1 Gertsenshtein channel at this nominal
parameter point). If any run diverges, **stop and debug** — this is a
solver regression.

## Sweep 2 — 2D (α, δₘ)

56 points (8 α × 7 δₘ), ~4 s total parallel. Submit to standard
`sapphire` queue:

```bash
scripts/hpc_shuttle.sh submit \
    --template scripts/hpc_templates/sapphire_cpu.sbatch \
    --name dp_2d \
    --nodes 1 \
    --time 01:00:00 \
    --cmd 'TIDAL_PARALLEL=112 bash examples/torsion_dark_photon/sweep_2d.sh'
```

## Sweep 3 — Monte Carlo (α, ξ, δₘ)

200-sample Latin-hypercube (`n-samples=200` in `sweep_mc.sh`). Same
submission shape as sweep 2.

```bash
scripts/hpc_shuttle.sh submit \
    --template scripts/hpc_templates/sapphire_cpu.sbatch \
    --name dp_mc \
    --nodes 1 \
    --time 01:00:00 \
    --cmd 'TIDAL_PARALLEL=112 bash examples/torsion_dark_photon/sweep_mc.sh'
```

## Monitoring and retrieval

```bash
scripts/hpc_shuttle.sh status           # one-shot, no polling
scripts/hpc_shuttle.sh tail <jobid>     # last 200 lines of slurm log
scripts/hpc_shuttle.sh pull <jobid>     # lightweight: CSV/JSON/figures only
```

## t_end-independence pairing (#238)

For publication-grade amplification claims, run each sweep at **two**
`t_end` values and join the CSVs via the `analyze_sweep.py --paired-dir`
option. Example for the 2D sweep:

```bash
# Primary sweep at t_end=50 (default in sweep_2d.sh)
scripts/hpc_shuttle.sh submit --name dp_2d_50 \
    --cmd 'TIDAL_PARALLEL=112 bash examples/torsion_dark_photon/sweep_2d.sh'
# Paired sweep at t_end=25 — override T_END in the script or
# invoke `tidal sweep` directly with --t-end 25 --output ...sweep_2d_t25
```

Then locally after pulling both:

```bash
python examples/torsion_dark_photon/analyze_sweep.py \
    examples/data/torsion_dark_photon_sweep_2d \
    --paired-dir examples/data/torsion_dark_photon_sweep_2d_t25
```

Any run with `P(t_end) / P(t_end/2) > 4.1` is flagged as a tachyonic
amplification artifact (`flagged_runs.csv` with reason `super-sin²`).

## Known cosmetic warnings

- `UserWarning: mass_matrix/coupling_matrix inconsistent with identity
  operator terms` — JSON loader warning, unrelated to #256.
- Modal solver `max Re(λ) ≈ 1` warning — spurious gauge/null modes
  that are correctly frozen downstream by the `|λ| > 1e12` gauge
  filter and `_suppress_tachyonic_noise`. See
  `docs/tex/troubleshooting.tex` §"Rank-Deficient Kinetic from
  Trace-Projection Lagrangians".

## References

- `SMOKE_TEST.md` — local smoke-sweep evidence for the post-#256 fix.
- `memory/reference_csd3_hpc.md` — CSD3 SSH / billing / partition
  authoritative facts.
- `memory/project_pgt_sweep_campaign.md` — running campaign log.
