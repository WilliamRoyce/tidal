"""``tidal sample`` — Bayesian inference via Monte Carlo and nested sampling.

Wraps the simulation + measurement pipeline as a likelihood function
for parameter estimation and model comparison.

Supports:
- Simple Monte Carlo: ``--method mc --n-samples 100``
- Nested sampling (PolyChord): ``--method nested --nlive 100``
- User-specifiable priors: ``--prior "alpha=uniform:0.01:10"``
- Hard constraints: ``--constraint "xi > 0"``
- Parallel evaluation: ``--parallel N``

References
----------
Handley, W. et al. (2015) "PolyChord", MNRAS 453(4).
Handley, W. (2019) "anesthetic", JOSS 4(37).
Skilling, J. (2004) "Nested Sampling", AIP Conference Proceedings 735.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.inference._results import InferenceResult


def sample_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0915
    """Entry point for ``tidal sample``."""
    from pathlib import Path

    from tidal.cli._console import error_with_hint

    # --- Validate required arguments ---
    json_path = getattr(args, "json_path", None)
    if not json_path:
        error_with_hint(
            "No JSON specification provided.",
            ["Usage: tidal sample <spec.json> --prior 'NAME=DIST:ARGS' ..."],
        )
        return 1

    spec_path = Path(json_path)
    if not spec_path.exists():
        error_with_hint(
            f"File not found: {spec_path}",
            ["Check the path to the JSON equation specification."],
        )
        return 1

    # --- Parse priors (per-parameter --prior + optional --joint-prior) ---
    prior_specs: list[str] = getattr(args, "prior", []) or []
    joint_prior_specs: list[str] = getattr(args, "joint_prior", []) or []
    if not prior_specs and not joint_prior_specs:
        error_with_hint(
            "No priors specified.",
            [
                "Use --prior 'NAME=DIST:LOW:HIGH' (repeatable) or",
                "--joint-prior 'names=...;type=cubed_sphere;M=K;face=F;"
                "sub=S;r_lo=L;r_hi=H' (cubed-sphere joint prior)",
                "Distributions for --prior: uniform, log_uniform, normal, "
                "arctan_uniform",
                "Example: --prior 'alpha=uniform:0.01:10'",
            ],
        )
        return 1

    from tidal.inference._prior import (
        Prior,
        RadialAngularPrior,
        parse_joint_prior,
        parse_prior,
        prior_param_names,
        total_prior_ndim,
        validate_prior_names,
    )

    if len(joint_prior_specs) > 1:
        error_with_hint(
            "Only one --joint-prior per invocation is supported in this release.",
            [
                "Multiple cubed-sphere joint priors (e.g. one per spin "
                "sector) are deferred — the supervisor's directive is one "
                "monolithic sphere over all BSM couplings."
            ],
        )
        return 1

    try:
        per_param_priors: list[Prior] = [parse_prior(s) for s in prior_specs]
        joint_priors: list[RadialAngularPrior] = [
            parse_joint_prior(s) for s in joint_prior_specs
        ]
    except ValueError as e:
        error_with_hint(
            str(e),
            [
                "Check --prior format: NAME=DIST:ARG1:ARG2",
                "Check --joint-prior format: 'names=N1,N2,...;"
                "type=cubed_sphere;M=K;face=F;sub=S1_S2_...;r_lo=L;r_hi=H'",
            ],
        )
        return 1

    priors: list[Prior | RadialAngularPrior] = [
        *joint_priors,
        *per_param_priors,
    ]

    try:
        validate_prior_names(priors)
    except ValueError as e:
        error_with_hint(
            str(e),
            ["--prior names must be disjoint from --joint-prior names"],
        )
        return 1

    # --- Parse constraints ---
    constraint_specs: list[str] = getattr(args, "constraint", []) or []
    from tidal.inference._constraints import ConstraintError, ConstraintSet

    try:
        constraints = ConstraintSet.from_strings(constraint_specs)
    except ConstraintError as e:
        error_with_hint(str(e), ["Check --constraint syntax: 'xi > 0'"])
        return 1

    # --- Parse likelihood ---
    likelihood_spec: str | None = getattr(args, "likelihood", None)
    if not likelihood_spec:
        error_with_hint(
            "No likelihood specified.",
            [
                "Use --likelihood 'METRIC:TYPE[:ARGS]'",
                "Types: maximize, gaussian:TARGET:SIGMA, threshold:MIN_VALUE",
                "Example: --likelihood 'P_max:maximize'",
            ],
        )
        return 1

    from tidal.inference._likelihood import parse_likelihood

    baseline_formula: str | None = getattr(args, "baseline_formula", None)
    # v3 architecture: --soft-floor-noise tunes the Normal(0, sigma) noise
    # on the soft penalty floor.  (There is no probe gate to configure:
    # rejection on tachyonic growth was abandoned as policy and the
    # --gated flag removed in v0.49.6 — see docs/tex/stability_probe.tex.)
    soft_floor_noise: float = float(getattr(args, "soft_floor_noise", 1.0))
    from tidal.inference._likelihood import SOFT_FLOOR_LOGL

    soft_floor_logl: float = float(getattr(args, "soft_floor_logl", SOFT_FLOOR_LOGL))
    # Read --seed here rather than at the "Common settings" block below: the
    # likelihood needs it so soft-floor noise is reproducible (issue #408).
    seed: int = getattr(args, "seed", 42)
    try:
        likelihood_config = parse_likelihood(
            likelihood_spec,
            baseline_formula=baseline_formula,
            soft_floor_noise_sigma=soft_floor_noise,
            soft_floor_logl=soft_floor_logl,
            noise_seed=seed,
        )
    except ValueError as e:
        error_with_hint(str(e), ["Check --likelihood format: METRIC:TYPE[:ARGS]"])
        return 1

    # --- Parse measurements ---
    measure_str: str | None = getattr(args, "measure", None)
    if measure_str:
        measurements = set(measure_str.split(","))
    else:
        # Auto-detect from likelihood metric
        metric = likelihood_config.metric
        if metric.startswith("P_") or metric in {"P_max", "P_final"}:
            measurements = {"conversion", "peak_conversion"}
        elif metric.startswith("E_") or metric == "max_energy_error":
            measurements = {"energy", "conservation"}
        elif metric == "L_mix":
            measurements = {"mixing"}
        else:
            measurements = {"summary"}

    source = tuple(args.source.split(",")) if getattr(args, "source", None) else None
    target = tuple(args.target.split(",")) if getattr(args, "target", None) else None
    threshold = getattr(args, "energy_threshold", 1e-3)

    # --- Validate output ---
    output_dir = getattr(args, "output", None)
    if not output_dir:
        error_with_hint(
            "No output directory specified.",
            ["Use --output DIR to specify where results are saved."],
        )
        return 1
    output_path = Path(output_dir)

    # Per-tile output redirect for cubed-sphere joint priors: chains land
    # in <output>/<face_label>_tile<sub_tile>/ so the atlas plotter can
    # pool tiles transparently.  Convention matches the supervisor's
    # psalter cell-directory layout.
    if joint_priors:
        from tidal.inference._sphere import cell_label

        jp = joint_priors[0]
        cell_dir = cell_label(jp.face_idx, jp.sub_tile)
        output_path /= cell_dir
        output_path.mkdir(parents=True, exist_ok=True)

    # --- Common settings ---
    method: str = getattr(args, "method", "mc")
    n_workers: int | None = getattr(args, "parallel", None)
    quiet: bool = getattr(args, "quiet", False)

    param_names = prior_param_names(priors)

    # --- Print summary ---
    if not quiet:
        print(f"=== tidal sample ({method}) ===")
        print(f"  Spec: {spec_path}")
        print(f"  Parameters: {', '.join(param_names)}")
        for p in priors:
            if isinstance(p, RadialAngularPrior):
                radial = (
                    f"fixed r={p.r_lo}"
                    if p.is_fixed_radius
                    else f"log_uniform(r={p.r_lo}..{p.r_hi})"
                )
                print(
                    f"    [{', '.join(p.names)}] ~ "
                    f"cubed_sphere(face={p.face_idx}, sub={p.sub_tile}, "
                    f"M={p.M}) x {radial}"
                )
            else:
                print(f"    {p.name} ~ {p.distribution}({p.low}, {p.high})")
        if constraints:
            print(f"  Constraints: {len(constraints)}")
            for expr in constraints.expressions:
                print(f"    {expr}")
        print(
            f"  Likelihood: {likelihood_config.metric} ({likelihood_config.likelihood_type})",
        )
        print(f"  Measurements: {', '.join(sorted(measurements))}")
        print(f"  Output: {output_path}")

    # --- Run inference ---
    if method == "mc":
        if joint_priors:
            error_with_hint(
                "--joint-prior is not supported with --method mc in this release.",
                [
                    "Use --method nested for cubed-sphere joint priors",
                    "(MC support for the joint sampler is deferred — the v3 "
                    "campaign use case is PolyChord nested sampling)",
                ],
            )
            return 1
        n_samples = getattr(args, "n_samples", 100)
        if not quiet:
            print(f"  Samples: {n_samples}")
            if n_workers:
                print(f"  Workers: {n_workers}")
            print()

        from tidal.inference._mc import run_monte_carlo

        # MC path is gated above to per-parameter priors only (joint
        # priors fail-fast); pass per_param_priors explicitly so the
        # type is `list[Prior]`.
        result = run_monte_carlo(
            priors=per_param_priors,
            likelihood_config=likelihood_config,
            base_args=args,
            spec_path=spec_path,
            measurements=measurements,
            source=source,
            target=target,
            threshold=threshold,
            n_samples=n_samples,
            constraints=constraints,
            n_workers=n_workers,
            seed=seed,
            temp_dir=output_path / "_runs",
            quiet=quiet,
        )

    elif method == "nested":
        nlive = getattr(args, "nlive", 100)

        # Auto-scale nlive if requested
        nlive_auto = getattr(args, "nlive_auto", None)
        if nlive_auto is not None:
            from tidal.inference._nested import recommend_nlive

            nlive = recommend_nlive(total_prior_ndim(priors), nlive_auto)
            if not quiet:
                print(f"  nlive auto ({nlive_auto}): {nlive}")

        if not quiet:
            print(f"  Live points: {nlive}")
            print("  Sampler: polychord")
            if n_workers:
                print(f"  Workers: {n_workers}")
            print()

        from tidal.inference._likelihood import SimulationLikelihood
        from tidal.inference._nested import run_nested_sampling
        from tidal.inference._prior import build_prior_transform

        likelihood_fn = SimulationLikelihood(
            base_args=args,
            spec_path=spec_path,
            param_names=param_names,
            measurements=measurements,
            source=source,
            target=target,
            threshold=threshold,
            likelihood_config=likelihood_config,
            temp_dir=output_path / "_runs",
        )

        # Collect optional PolyChord-specific settings (pass-through kwargs)
        ns_kwargs: dict[str, object] = {
            "output_dir": str(output_path / "_chains"),
        }
        num_repeats = getattr(args, "num_repeats", None)
        if num_repeats is not None:
            ns_kwargs["num_repeats"] = num_repeats
        precision_criterion = getattr(args, "precision_criterion", None)
        if precision_criterion is not None:
            ns_kwargs["precision_criterion"] = precision_criterion
        if getattr(args, "no_clustering", False):
            ns_kwargs["do_clustering"] = False
        if getattr(args, "read_resume", False):
            ns_kwargs["read_resume"] = True
        max_ndead = getattr(args, "max_ndead", None)
        if max_ndead is not None:
            ns_kwargs["max_ndead"] = max_ndead

        result = run_nested_sampling(
            log_likelihood=likelihood_fn,
            prior_transform=build_prior_transform(priors),
            ndim=total_prior_ndim(priors),
            param_names=param_names,
            sampler="polychord",
            nlive=nlive,
            n_workers=n_workers,
            quiet=quiet,
            **ns_kwargs,
        )
        # Record priors in metadata so post-hoc analysis (e.g. correct
        # marginal D_KL per parameter, #308) can transform each column
        # into the space where its prior is uniform.  Mixed priors are
        # serialized by type so the loader can reconstruct each component.
        from tidal.inference._prior import to_record

        result.metadata["priors"] = [to_record(p) for p in priors]
        # Persist likelihood mode so the corner-plot title can label the
        # extremal logL correctly: "max A" for maximize, "max suppression"
        # for minimize, etc.  Without this the same number means different
        # things and confuses readers (suppress's high logL maps to
        # tiny P_max, not large amplification).
        result.metadata["likelihood_type"] = likelihood_config.likelihood_type
        result.metadata["likelihood_metric"] = likelihood_config.metric
        # Persist enough of the simulation context to recompute the
        # baseline formula (and therefore amplification A = P_max / P_GR)
        # post-hoc — needed for the corner-plot's derived A column after
        # results are pulled from HPC where the original CLI args aren't
        # available.  See tidal/inference/_visualize.py.
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )
        from tidal.symbolic import load_equation_system as _load_spec_for_meta

        try:
            spec_for_params = _load_spec_for_meta(spec_path)
            persisted_params = _parse_params(
                list(getattr(args, "param", []) or []), spec_for_params
            )
        except Exception:  # noqa: BLE001
            persisted_params = {}
        sim_params: dict[str, object] = {
            "t_end": float(getattr(args, "t_end", 0.0) or 0.0),
            "fixed_params": persisted_params,
            "baseline_formula": likelihood_config.baseline_formula,
        }
        result.metadata["simulation_params"] = sim_params

        # Post-hoc prior stability sweep: PolyChord drops -inf samples
        # before they enter the chain, so the unstable region is invisible
        # in the chain.  Sample the prior independently and run only the
        # cheap eigenvalue check (~1 ms per draw, no simulation) to
        # generate a side file for corner-plot overlay.  See
        # tidal/inference/_prior_stability.py for the rationale.
        if source and target:
            try:
                from tidal.inference._prior_stability import (
                    run_prior_stability_sweep,
                )

                rej_path = output_path / "_rejected_prior.csv"
                run_prior_stability_sweep(
                    base_args=args,
                    spec_path=spec_path,
                    param_names=param_names,
                    prior_transform=build_prior_transform(priors),
                    source=source,
                    target=target,
                    output_path=rej_path,
                    n_samples=getattr(args, "prior_sweep_samples", 5000),
                    seed=seed,
                    quiet=quiet,
                )
                result.metadata["rejected_prior_path"] = str(rej_path)
            except Exception as exc:  # noqa: BLE001
                if not quiet:
                    print(f"  Prior stability sweep skipped: {exc}")

    else:
        error_with_hint(
            f"Unknown method '{method}'.",
            ["Use --method mc or --method nested"],
        )
        return 1

    # --- Save results (rank 0 only — all MPI ranks reach here) ---
    # On HPC login nodes mpi4py imports fine but COMM_WORLD initialization
    # aborts (PMI2 not available outside a SLURM allocation).  Catch any
    # error and treat as rank 0 — the only consequence is that all ranks
    # in a real MPI run would write the file (rank-0 logic is a wallclock
    # optimization, not correctness-critical).
    try:
        from mpi4py import MPI  # type: ignore[import-untyped]

        mpi_rank: int = int(MPI.COMM_WORLD.Get_rank())  # type: ignore[reportUnknownArgumentType]
    except Exception:  # noqa: BLE001
        mpi_rank = 0
    if mpi_rank != 0:
        return 0

    result.save(output_path)
    if not quiet:
        print()
        _print_summary(result)
        print(f"\nResults saved to: {output_path}")

    # --- Optional plots ---
    if getattr(args, "corner", False):
        from tidal.inference._visualize import plot_corner

        corner_path = output_path / "corner.png"
        plot_corner(result, corner_path)
        if not quiet:
            print(f"Corner plot: {corner_path}")

    if getattr(args, "trace", False):
        from tidal.inference._visualize import plot_trace

        trace_path = output_path / "trace.png"
        plot_trace(result, trace_path)
        if not quiet:
            print(f"Trace plot: {trace_path}")

    if getattr(args, "importance", False):
        from tidal.inference._importance import format_importance_table
        from tidal.inference._visualize import plot_importance

        importance_path = output_path / "importance.png"
        try:
            imp = result.parameter_importance()
            plot_importance(imp, importance_path)
            if not quiet:
                # Print the annotated table alongside the chart: the chart
                # alone carried none of the consistency warnings (floor,
                # fallback, superadditivity) — the #433 review gap where
                # `--importance` was the one surface with zero caveats.
                print(format_importance_table(imp))
                print(f"Importance plot: {importance_path}")
        except (ImportError, ValueError) as e:
            if not quiet:
                print(f"  (importance plot skipped: {e})")

    # --- Optional analysis ---
    if getattr(args, "analyze", False) and result.method == "nested":
        try:
            from tidal.inference._importance import format_importance_table

            imp = result.parameter_importance()
            if not quiet:
                print(format_importance_table(imp))
        except (ImportError, ValueError) as e:
            if not quiet:
                print(f"  (importance analysis skipped: {e})")

    return 0


def _print_summary(result: InferenceResult) -> None:
    """Print inference summary statistics."""
    r = result

    print(f"  Method: {r.method}")
    print(f"  Samples: {r.n_samples}")
    print(f"  ESS: {r.effective_sample_size():.0f}")

    if r.log_evidence is not None:
        print(f"  log Z = {r.log_evidence:.2f} +/- {r.log_evidence_err:.2f}")

    print("  MAP estimate:")
    for name, val in r.best().items():
        print(f"    {name} = {val:.6g}")

    print("  Posterior mean:")
    for name, val in r.posterior_mean().items():
        print(f"    {name} = {val:.6g}")

    ci = r.credible_interval(0.95)
    print("  95% credible intervals:")
    for name, (lo, hi) in ci.items():
        print(f"    {name}: [{lo:.6g}, {hi:.6g}]")
