# V3 Phase C: cubed-sphere coupling-space chart — reference

This document records the implementation of the v3 plan's **Phase C — sphere-cube radial-angular prior**, which the original v3 plan deferred pending a supervisor reply. The supervisor's package `psalter` (sent 10 May 2026; gitignored, treated as opaque local reference only) and his subsequent email replies pinned the conventions; the implementation landed in [`8f95157`](https://github.com/.../commit/8f95157) and [`99da4b2`](https://github.com/.../commit/99da4b2). This doc serves as the durable handoff so subsequent work on the v3 plan can drop the now-obsolete items.

## What landed

| Module | Purpose |
| --- | --- |
| [tidal/inference/_sphere.py](../tidal/inference/_sphere.py) | Cubed-sphere geometry: face indexing, sub-tile bounds, gnomonic projection, optional rotation. ~200 LOC, pure numpy. |
| [tidal/inference/_prior.py](../tidal/inference/_prior.py) `RadialAngularPrior` | Joint prior over `(r, theta_hat)`: magnitude × cubed-sphere direction. Mixed with per-coupling `Prior` via the same `build_prior_transform`. |
| [tidal/inference/_atlas.py](../tidal/inference/_atlas.py) | Atlas plot: pools per-tile chains under `<survey>/<face_label>_tile<sub>/` and renders a 2N-face panel grid. |
| [tidal/cli/_sample.py](../tidal/cli/_sample.py) `--joint-prior` | CLI flag for cubed-sphere joint priors; per-tile output dir; mutually exclusive with `--method mc`; mixable with `--prior`. |
| [tidal/cli/_plot_command.py](../tidal/cli/_plot_command.py) `--type atlas` | CLI dispatch for atlas plot. |
| [tests/test_sphere_geometry.py](../tests/test_sphere_geometry.py), [tests/test_radial_angular_prior.py](../tests/test_radial_angular_prior.py), [tests/test_atlas_plot.py](../tests/test_atlas_plot.py) | 102 new tests covering geometry, prior, CLI parsing, atlas smoke. |

## Geometry primer (load this mental model first)

The cubed-sphere chart is initially confusing because the face panels' axes are not individual coupling values. This 6-step walkthrough loads the picture; the rest of this doc and the code assume it.

**1. What the "sphere" is.** The N couplings form a vector `c = (c_1, ..., c_N) ∈ R^N`. We re-express it as `(magnitude, direction)`: `r = sqrt(c_1² + ... + c_N²)` is the length, `theta_hat = c/r` is the unit-length direction. Every direction `theta_hat` is one point on the unit sphere `S^(N-1)` — a curved (N-1)-dim surface in N-dim Euclidean space (think globe surface for N=3). `r` and `theta_hat` are independent: same `theta_hat` but bigger `r` = same *relative* coupling pattern at higher overall strength; same `r` but different `theta_hat` = same overall strength but different relative pattern. Per the supervisor's reply, `N` here is the count of *non-standard / BSM* couplings only; EH and Maxwell prefactors live outside the sphere (more below).

**2. The cube wrapped around the origin.** Imagine a cube `[-1, 1]^N` centered at the origin. Its surface has `2N` faces (12 for N=6). Each face is the (N-1)-dim slab where one axis is pinned at `±1` and the others vary in `[-1, 1]`. Face 1 ↑ is `(+1, u_2, u_3, ..., u_N)`; Face 1 ↓ is `(-1, u_2, ...)`; Face 2 ↑ is `(u_1, +1, u_3, ..., u_N)`; and so on through Face N ↓.

**3. Radial projection cube → sphere (gnomonic).** For each point on the cube surface, draw a ray from the origin through it and record where it hits the unit sphere. That hit is the corresponding point on `S^(N-1)`. Each face's `[-1, 1]^(N-1)` cube patch becomes a quasi-rectilinear chart of the sphere region "where coupling k dominates with sign ±". The 2N face charts together cover the whole sphere with mostly-uniform spacing — that's the cubed-sphere's selling point over latitude/longitude (which crowds at the poles).

**4. What each face panel shows.** For Face 1 ↑, the panel's axes are the face-local cube coordinates `chi_1^{1↑}` through `chi_{N-1}^{1↑}`, each in `[-1, 1]`. A 2D KDE filled contour `chi_2^{1↑}` vs `chi_1^{1↑}` shows the marginal posterior at positions on the cube face — i.e., what direction-on-the-sphere-near-the-+x_1-pole the sampler concluded was favoured. Reading a contour peak at `(chi_1, chi_2, chi_3, chi_4, chi_5) = (0.3, -0.5, 0.1, 0.2, 0)` on Face 1 ↑ means: build the cube vector `(+1, 0.3, -0.5, 0.1, 0.2, 0)`, normalise to unit length, that's the favoured `theta_hat` — corresponding to physical-coupling *relative pattern* `(0.847, 0.254, -0.423, 0.085, 0.169, 0)`.

**5. What the axes' "scales" mean.** The face-local axes are NOT individual coupling values — they are cube chart coordinates with fixed range `[-1, 1]` per axis. To recover individual physical couplings `c_i` from a sample at position `chi` on Face k ↑/↓:

1. Build the cube vector `v` with `v_k = ±1` and `v_j = chi_j^{k±}` for the other N-1 slots.
2. Normalise: `theta_hat = v / |v|`.
3. Multiply by the magnitude: `c = r · theta_hat`.

So `c_i = r · v_i / |v|`. The chi axes encode the *relative* sign-and-weight of each coupling; the dominant coupling's contribution to `theta_hat` is `±1/|v|`. The inverse direction (chain → face-local chi) is implemented in `tidal.inference._atlas._physical_to_face_local`.

**6. What varying `r` does.** `r` is sampled independently (a separate 1D dimension of the joint prior, e.g. `log_uniform(1e-3, 1e3)`). Varying `r` rescales *all* physical couplings together by the same factor. Every face panel still shows positions on the *unit* sphere; those positions don't change with `r`. `r` is plotted separately — typically as a 1D histogram of the posterior over `log r`. So a full posterior under this scheme is: `2N` face panels (showing the direction posterior over `theta_hat`) + one 1D `log r` histogram (showing the magnitude posterior). Face panels capture the *relative coupling pattern*; the `log r` histogram captures the *overall coupling strength*.

**Concrete example for N=6.** Suppose the chain finds two distinct regions:

- Peak A: `c = (1.0, 0.2, -0.3, 0, 0, 0)` → `r_A ≈ 1.063`, `theta_hat_A ≈ (0.941, 0.188, -0.282, 0, 0, 0)`. Dominant axis +x_1 → lives on Face 1 ↑ at face coords `(0.2, -0.3, 0, 0, 0)`. KDE peak there.
- Peak B: `c = (0, 0, 0, 0, 0, 5.0)` → `r_B = 5.0`, `theta_hat_B = (0, 0, 0, 0, 0, +1)`. Dominant axis +x_6 → lives on Face 6 ↑ at face coords `(0, 0, 0, 0, 0)` (face centre). KDE peak at the centre.
- The `log r` posterior is bimodal: hump at `log 1.063 ≈ 0` and hump at `log 5 ≈ 0.7`.

You read the science off the joint plot: **which face panels show density** = which coupling-direction combinations the sampler likes; **the 1D `log r` plot** = at what overall strength. Faces with no density = sphere regions the chain ruled out.

## Convention map (matches `psalter.tile`)

The supervisor's `psalter` package — sent as a tarball, gitignored, *not* vendored — pins the exact conventions. We re-implemented from a clean reading of standard cubed-sphere geometry (Ronchi et al. 1996); psalter is the reference for naming so per-tile outputs are interchangeable for cross-checking.

- **Face indexing**: 1-indexed `[1, 2N]`. `face 1 = +x_1`, `face 2 = -x_1`, `face 3 = +x_2`, .... `face_label` returns `"01p"`, `"03m"`; `face_label_math` returns `r"1\uparrow"`, `r"6\downarrow"`.
- **Cube → sphere map**: gnomonic — `theta = (axis_vec(face) + embed @ u) / |·|`, then optional `Q.T @ x` pre-rotation.
- **Sub-tiles**: per-face axes M-subdivided into M^(N-1) cells indexed by 1-indexed tuple `(s_1, ..., s_{N-1})`; bounds `u_lo = -1 + 2(s-1)/M`, `u_hi = -1 + 2s/M`.
- **Cell directory layout**: `<face_label>_tile<s_1>_<s_2>_...` — e.g. `01p_tile2_3` is face 1+, sub-tile (2, 3). Built via `tidal.inference._sphere.cell_label`.
- **Q rotation**: identity by default. `random_rotation(N, seed)` is available (proper rotation via QR of N(0, I), det = +1) for the "trick for more readable plots" use case the supervisor described — not load-bearing physics.

## What goes on the sphere vs off

Per the supervisor's two-message reply (10 May 2026) — initial directive then explicit retraction:

- **Standard kinetic prefactors (Einstein-Hilbert, Maxwell)** are **carved out** of the sphere. This is a deliberate dimensionality reduction (walltime win); the physical roles of these couplings are well-established, so spending sphere dimensions on them costs nothing scientific.
- **All non-standard / BSM operators** go on the cubed-sphere on equal footing. No further grouping by symmetry class, SPO sector, or parity. Physical meaning of each BSM coupling emerges from the sampling outcomes.

In current TIDAL practice, the carve-out is *already realised*: `kappa` lives in the `[constants]` block of theory TOMLs fixed at 1.0, and the Maxwell term is hard-coded as `-1/4 F_{ab} F^{ab}` with no free prefactor. Current campaign scripts (e.g. `scripts/hpc_submit_drafts/d1_*.sh`) only `--prior` over the BSM couplings. The cubed-sphere prior shipped here therefore operates over *exactly* the set of couplings that already get `--prior` flags.

## Drop-in replacement for current `--prior` blocks

Existing d1 campaign uses:

```text
--prior "alpha1=uniform:-1:1" \
--prior "alpha2=uniform:-2:2" \
--prior "alpha3=log_uniform:0.05:2" \
--prior "delta1=uniform:-2:2"
```

The cubed-sphere drop-in is one flag:

```text
--joint-prior "names=alpha1,alpha2,alpha3,delta1;type=cubed_sphere;M=2;face=1;sub=1_1_1;r_lo=1e-3;r_hi=1e3"
```

— a 4-coupling sphere with `r ~ log_uniform(1e-3, 1e3)`, `theta_hat` on a 2-subdivided cubed-sphere of `S^3` (N=4 → `2N=8` faces, M^(N-1) = 8 sub-tiles per face, 64 cells total). For one campaign, the survey driver fans out 64 `tidal sample --joint-prior ...` invocations (one per cell) via `scripts/hpc_shuttle.sh` array jobs.

Per-tile output goes to `<output>/<face_label>_tile<sub>/...`. `tidal plot <output> --type atlas` then pools all tiles transparently.

## v3 plan items the parallel agent can drop

The v3 plan at `/home/vscode/.claude/plans/binary-snacking-chipmunk.md` lists Phase C as DEFERRED and earmarks specific tasks. With this implementation landed, those tasks are now closed:

1. **Phase C entirely.** `tidal.inference._sphere`, `RadialAngularPrior`, `tidal sample --joint-prior`, `tidal plot --type atlas` are live. Drop the deferred Phase C from active tracking.
2. **Phase C.2 supervisor email** — closed by his 10 May reply (saved at `docs/meetings/2026-05-10_supervisor_reply_grouping.md`). No follow-up needed.
3. **Phase A.0 `docs/V3_PHASE_C_DESIGN.md`** — superseded by *this* doc. Either delete the placeholder task or repoint it to here.

Items the parallel agent should **keep** on their plan:

- Phase A.2 per-param `arctan_uniform` rewrites — separate stopgap; both compactifications coexist, campaigns choose `--prior` (A.2 path) or `--joint-prior` (this) per run.
- Phase A.3 corner-plot upper-triangle hide at `tidal/inference/_visualize.py:758-789` — orthogonal to atlas plot; still relevant for non-cubed-sphere posteriors.
- Phases A.0 GitHub issues, A.1 soft-penalty refactor, A.4 tests, A.5 D1 v1 replay, A.6 commit/version bump, A-γ γ_conversion, B campaign re-runs, D manuscript, E localised geometry — untouched.

## Verification

| Step | Command | Status |
| --- | --- | --- |
| Geometry tests (48) | `uv run pytest tests/test_sphere_geometry.py -x -q` | PASS |
| Prior tests (29) | `uv run pytest tests/test_radial_angular_prior.py -x -q` | PASS |
| Atlas tests (25) | `uv run pytest tests/test_atlas_plot.py -x -q` | PASS |
| Full suite (2200) | `uv run pytest tests/ -x -q` | PASS |
| Type / lint / format | `uv run pyright tidal/inference/_sphere.py tidal/inference/_prior.py tidal/inference/_atlas.py tidal/cli/_sample.py tidal/cli/_plot_command.py && uv run ruff check && uv run ruff format` | PASS |

## Reference design (psalter, not vendored)

The supervisor's package is read only via `tar -xzOf psalter.tar.gz <member>` — never extracted to the working tree. The relevant reference modules are `psalter/_tile/geometry.py` and `psalter/tile.py` (geometry conventions); `psalter/_sample/sample.py` (per-tile NS driver pattern); `psalter/_plot/core.py` (atlas renderer pattern).

Our implementation is independent code and not a derivative of psalter; the value of psalter as a reference is that it pins the *naming and orientation conventions* — so a TIDAL chain in `01p_tile2_3/` and a psalter chain in the same directory layout describe the same sphere region, simplifying any future cross-validation. We do not redistribute psalter; neither does this codebase.
