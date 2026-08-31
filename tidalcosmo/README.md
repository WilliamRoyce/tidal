# `tidalcosmo/` — the Cobaya extension (scaffold)

> **Preliminary — planning-stage scaffold (H4, 2026-08-31).** These directories and their stated
> responsibilities come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. They record why each boundary was drawn so a
> later reader can weigh it — **they are expected to be revised or replaced.** Changing one does
> not require re-litigating H4.

**There is no code here yet.** This tree is directories and READMEs only. Deliberately no
`__init__.py` and no Python: an empty importable package would enter pyright, ruff and coverage
with nothing in it.

**Full design:** [`docs/cosmology/repo_reshape.md`](../docs/cosmology/repo_reshape.md).
**Program record:** [`docs/COSMOLOGY_PROGRAM.md`](../docs/COSMOLOGY_PROGRAM.md).
**New to the field?** [`docs/cosmology/primer.md`](../docs/cosmology/primer.md), then §0 of the
design document, which is a glossary written for a theoretical physicist new to computational
cosmology.

## What this package is

A **Cobaya extension**: a `Theory` class that chains off upstream CAMB, evolves a candidate
Lagrangian's perturbations as spectators on a ΛCDM background, and returns CMB observables.

## Three rules that bind every session

1. **`tidalcosmo` is a placeholder for `tidal`.** The end state is: build this alongside the
   legacy `tidal/` tree, migrate capability by capability, delete `tidal/`, then
   `git mv tidalcosmo tidal`. So: **no cosmology-specific naming inside this tree** (every
   directory below renames cleanly), assume the final console script is `tidal` and the final
   dotted path is `tidal.SpectatorTheory`, and **do not publish or circulate anything naming
   `tidalcosmo` outside the project** — an external `theory: {tidalcosmo.SpectatorTheory: …}`
   breaks on the day we rename.

2. **New code never imports old code.** No adapters, no shims, no subprocess calls into the
   legacy CLI. This is enforced by a hygiene test over `tidalcosmo/` and `tests_cosmo/`. The
   legacy tree is a **test oracle with an expiry date**, so it is captured as *committed golden
   data* under `tests_cosmo/data/oracles/`, never as a live dependency — the moment a test
   imports or shells out to legacy, legacy becomes undeletable.

3. **Conform to external conventions natively.** This package adopts the naming, formats, gauge
   conventions and interchange of the tools it interoperates with — **CAMB and PSALTer** — from
   the start. Legacy TIDAL notation has no claim here: there is no backward compatibility to
   preserve, and no conversion layer to build. See design §2.8. One direct consequence: the
   `derive` gate is *semantic equivalence* to frozen legacy specs, **not** a byte diff.

## Layout, and why it is two tiers

**The surface is organized by what a user names in their Cobaya YAML.** Cobaya requires a
component's defaults to live in `<ClassName>.yaml` *beside the module defining the class*, so the
directory is addressed by the config file rather than freely chosen.

| Directory | Responsibility |
| --- | --- |
| `spectator/` | the `SpectatorTheory` Cobaya component and its per-channel YAML presets |
| `background/` | the CAMB seam; `protocol.py` is H3's backend swap point. Seam A (tabulated `(a, ρ, P)`) is optional, off by default, and blocked on the CAMB fork (#498) |
| `likelihoods/` | only what the ecosystem does not already supply |
| `presets/` | runnable ladder configurations, O0 onward |

**The engine is organized by stage of the calculation**, named after DISCO-EB so a cosmologist
already knows what is in each file.

| Directory | Responsibility |
| --- | --- |
| `config/` | frozen dataclass settings bundles (PSALTer's idiom) — the fix for legacy's `argparse.Namespace`-as-config |
| `spec/` | the equation-spec interchange contract, carrying the declared gauge as metadata |
| `derive/` | the Wolfram front end — a port **and** a substantial FRW extension |
| `coefficients/` | symbolic coefficient → numeric callable of `(η, k)` |
| `perturbations/` | assembly of `M(η, k)`: our fields **plus the standard modes they couple to** |
| `solver/` | WS3; two front-ends over a shared core. Internals are H3's |
| `observables/` | transfer functions, line-of-sight projection, rotations |
| `validity/` | honest per-run flags — `ΔN_eff`, amplitudes, growth impact, background-EOM residual |
| `spectrum/` | WS6 polology — a **fast gate on the sampling path**, not a corner |
| `diagnostics/` | post-processing over Cobaya/anesthetic output |
| `cli/` | thin: parse → config → library call. No physics |

**The one structural inversion relative to legacy:** the forward model is library code taking a
typed config, and the CLI and the Cobaya component are two thin callers of the same entry point.
In the legacy package the CLI *is* the config object, which is why four modules outside
`tidal/cli/` import private names back out of it.
