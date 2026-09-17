# `perturbations/` — research artifacts behind the R-C perturbation-tooling memo

<!-- cspell:words xPand xMAG xBrauer TraceFree Pitrou Umeh Helpin Bertschinger xPert normu Weitzenböck Gorji Valcarcel -->

**Status: lane in progress (dispatched 2026-09-15; started 2026-09-16).** Results are written up in
`docs/cosmology/perturbation_tooling.md` (R-C, #567); the planning trail is archived verbatim in
`docs/cosmology/rc_planning_record.md`. Nothing here runs in the pipeline. **Research code, not
production code.**

## Contents

| path | what it is | provenance |
| --- | --- | --- |
| `install_xpand.sh` → **promoted to `scripts/install-xpand.sh`** (2026-09-17, D-C) | additive install of xPand 0.4.4 into `Applications/xAct/xPand/` (sha256-pinned tarball, junk stripped, `INSTALLED_VERSION` stamp, additivity asserted) | ours |
| `install_xmag.sh` | additive install of xMAG and its chain (SymmetricFunctions, BrauerAlgebra, xBrauer, TraceFree) at pinned commits, `INSTALLED_COMMIT` stamps | ours |
| `manifest.sh` | two-layer sha256 manifest of the whole userbase (must-be-identical layer / allow-listed caches / the new directories) | ours |
| `run_lane.sh` | serial lane runner: kernel guard, offscreen front end, throwaway cwd, hard timeout, scrubbed transcript, sentinel verdict | ours |
| `wolfram/RCSetup.wl` | shared setup: loads xPand first, FRW slicing with selectable normal norm, verdict / normalization / digest / export helpers, `RCTry` (a Throw becomes `$Failed` with its message visible), `RCSort` (xPand's own commutation of the flat slice derivatives, needed after `VarD`) | ours |
| `wolfram/probe_load.wls` | load on the certified bundle; the paper's minimal example (arXiv:1302.6174 TeX :1350-1373) reproduced; label layout and gauge rules read back | ours; targets from Pitrou–Roy–Umeh |
| `wolfram/repro_a2_tensor_eom.wls` | A2 control: FRW tensor equation from the perturbed Einstein equations, mixed and lowered indices, background residual printed | ours |
| `wolfram/repro_a1_tensor_action.wls` | A1 headline: action → tadpole → second order → Euler–Lagrange after the split | ours |
| `wolfram/repro_b_mb_scalars.wls` | Ma & Bertschinger scalar equations in the conformal-Newtonian and synchronous gauges | ours; targets from `literature/astro-ph_9506072/9506072.tex` |
| `wolfram/probe_c_signature.wls` | C: xPand under the project's (+,−,−,−) signature — `SetSlicing` with `normu = +1`, the paper's minimal example against the signature map, the `ExtractComponents` time-projection sign, the fluid split with the right norm | ours; hypotheses from `xPand.m:1715, 3014, 3042` |
| `wolfram/probe_d_torsion.wls` | D: a rank-3 torsion perturbation beside `dg` — SVT rules transcribed from arXiv:2310.16007 (TeX :768-842) and 1804.09215 (:819-852), projections, the action route for the torsion sector, the derivative split, the torsionful `DefCovD` connection and legacy's `ChangeCurvature` rewrite; optional key `d5=1` for a homogeneous background mode | ours; parametrization from Aoki–Bahamonde–Gigante Valcarcel–Gorji |
| `wolfram/probe_d_limits.wls` | D, fresh kernel: the rank map of xPand's split for a rank-3 field (scalars, rank 1, rank 2, rank 3, rank 4) with a pure-metric sanity split between attempts; the torsionful `DefCovD` connection and legacy's `ChangeCurvature` route in the same clean sequence | ours |
| `wolfram/probe_e_xmag.wls` | E: xMAG 0.1.0 beside xPand — load, the Riemann–Cartan `DefCovD`, contortion vs legacy's identity, `ToDistortion` + `BreakDistortion` vs legacy's `ChangeCurvature` route, `DefConnectionPerturbation` beside `DefMetricPerturbation`, end-to-end through xPand with Probe D's rules, `StartInducedDecomposition` | ours; calling forms from xMAG's `Documentation/DefCovD_xMAG.nb` |
| `wolfram/probe_e2_xmag_induced.wls` | E (f) in isolation: xMAG's `StartInducedDecomposition` on a kernel with xPand loaded but no xPand slicing, on a general and on a Riemann–Cartan connection | ours |
| `wolfram/RCSetupCore.wl` | the package-independent half of the harness: verdicts, digests, WXF export, `RCTry` (a Throw becomes `$Failed`, messages stay visible), `RCNewMessages`, `RCStackOn`, `RCSigns`/`RCSetSigns`. Read **after** the package a step is about, which its own guard enforces | ours |
| `wolfram/tier1_xmag.wls` | xMAG loaded first and called as its author documents: his `DefCovD_xMAG.nb` cells replayed against his stored outputs, the curvature-sign experiment with controls, and the three calls R-C got wrong | ours; targets are the author's stored outputs |
| `wolfram/f2_hazards.wls` | the two "xPand breaks" cases separated (chart derivative vs a `Master`-less connection), the `Master -> g` twin control, the caller-side `$CovDs` wrapper with a wrong-guard control, the xCoba count, and the shared-kernel control | ours |
| `wolfram/f3_changecurvature.wls` | `ChangeCurvature` on a `Master`-less torsion connection: plain xTensor (control) vs xMAG loaded, with the inert `TorsionToDistortion` witness | ours |
| `wolfram/f5_sign_audit.wls`, `wolfram/f6_psalter_signs.wls` | the five xTensor sign globals in each state: plain, xPand, after `SetSlicing`, xMAG, PSALTer — the lines `conventions.md` quotes | ours |
| `mb_camb_symbolic.py`, `camb_symbolic_newtonian_2.0.4.txt` | CAMB's own scalar equations printed from `camb.symbolic` (camb 2.0.4) — the machine check of the MB → our-convention transcription | ours; equations are CAMB's |

## How to reproduce

All Wolfram steps are strictly serial (one kernel, machine-wide). Outputs go to gitignored
`third_party/perturbations_runs/`. Take the userbase manifest before the gates and after them.

```bash
bash scripts/research/perturbations/manifest.sh before
bash scripts/psalter/ensure_registered.sh            # registers, then verify --require-psalter; exit 0 required
bash scripts/install-xpand.sh                        # no kernel; no-op if 0.4.4 is present
bash scripts/research/perturbations/run_lane.sh probe_load
bash scripts/research/perturbations/run_lane.sh a2
bash scripts/research/perturbations/run_lane.sh b
bash scripts/research/perturbations/run_lane.sh a1
bash scripts/research/perturbations/run_lane.sh c
bash scripts/research/perturbations/run_lane.sh d
bash scripts/research/perturbations/run_lane.sh dl
bash scripts/research/perturbations/install_xmag.sh  # no kernel; refuses if any target exists
bash scripts/research/perturbations/run_lane.sh e
bash scripts/research/perturbations/run_lane.sh e2
bash scripts/research/perturbations/run_lane.sh t1   # xMAG as documented, vs the author's outputs
bash scripts/research/perturbations/run_lane.sh f2   # the two session hazards, with controls
bash scripts/research/perturbations/run_lane.sh f3   # ChangeCurvature on a Master-less connection
bash scripts/research/perturbations/run_lane.sh f5   # sign globals: plain, xPand, xMAG
bash scripts/research/perturbations/run_lane.sh f6   # sign globals: PSALTer (own kernel)
bash scripts/verify-wolfram-setup.sh --require-psalter
bash scripts/research/perturbations/manifest.sh after && bash scripts/research/perturbations/manifest.sh diff
uv run python scripts/research/perturbations/mb_camb_symbolic.py scripts/research/perturbations/camb_symbolic_newtonian_2.0.4.txt
```

## Working with a third-party xAct package — the protocol

Three of this lane's "the package is broken" findings were its own calling forms, and one was
a convention the package sets silently (memo §3.3, §8; the PSALTer #543 precedent is the same
pattern). Before any claim about a third-party package, in this order:

1. **Read its documentation where it actually lives.** The usage strings — `?Symbol` in a
   kernel, and the `::usage` block at the top of the `.m` file (`xPand.m:127-668`,
   `xMAG.m:120-341`), plus xTensor's usage for every validator the package forwards to. The
   notebooks with stored outputs (`xMAG/Documentation/*.nb`,
   `xPand/Documentation/English/xPandDoc.nb`, `xPand/Examples/*.nb`). The README, the
   `.History` file, the paper. Then quote, in a comment above each call, the usage line the
   call follows.
2. **Replay one of the author's own examples first** in a clean kernel, loaded in the
   author's order, and compare with his stored outputs (`tier1_xmag.wls` does this with
   `Get` on the notebook and `ToExpression` on the cell boxes). A mismatch there is a
   bundle-compatibility finding; a mismatch only in your own call is your call.
3. **Check the options the package's own text calls important** (for xMAG, `Master`) and
   every validator it forwards to xTensor (a derivative's postfix symbol must be exactly one
   character, `xTensor.m:6591`).
4. **Print the sign globals** with `RCSigns` at kernel start, after every `Needs` and after
   anything that might change them; never compare expressions computed in two kernels
   without both sign headers. `RCSetSigns` is the only way to set them (short names before
   xTensor loads create shadowing `Global` symbols).
5. **Never `$MessageList = {}`** — it is protected, and the reset silently does nothing. Use
   `RCNewMessages[tag, expr]`, which prints only the messages that call produced.
6. **A message that names your own argument is your bug** until proven otherwise, and a
   benign message is not a failure (`DefMetric::old` appears in the author's own successful
   run, which is why `RCTry` no longer uses `Check`).
7. **Reduce to the author's names and the smallest session** before writing anything up, and
   read the package's known issues first (`COSMOLOGY_PROGRAM.md:591`).

## Rules this directory follows

- **Route, not payload.** Upstream packages are fetched by pinned sha256 or commit into the
  userbase; nothing is vendored into the repository.
- **One Wolfram kernel at a time.** Every launch goes through `run_lane.sh` (guard, offscreen,
  throwaway cwd, timeout); a step is judged by its `RC_<STEP>_DONE` sentinel and the verdict
  lines it prints, never by its exit status. The lane hook cannot see a bare `.wls`, so the
  runner and both installers guard themselves.
- **Verdicts, not checksums.** `identical (SameQ)` or `proved-equal` (the canonical difference
  is zero) count as reproduced; `proved-different` and `could-not-decide` do not. Every
  reproduction carries a negative control that must come out `proved-different`. Digests are
  fingerprints only.
- **No machine-specific paths** in anything committed; transcripts are scrubbed at print time
  (`tests/test_repo_hygiene.py` scans tracked files).
