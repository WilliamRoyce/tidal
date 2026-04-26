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


def sample_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0914, PLR0915
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

    # --- Parse priors ---
    prior_specs: list[str] = getattr(args, "prior", []) or []
    if not prior_specs:
        error_with_hint(
            "No priors specified.",
            [
                "Use --prior 'NAME=DIST:LOW:HIGH' (repeatable)",
                "Distributions: uniform, log_uniform, normal, arctan_uniform",
                "Example: --prior 'alpha=uniform:0.01:10'",
            ],
        )
        return 1

    from tidal.inference._prior import parse_prior

    try:
        priors = [parse_prior(s) for s in prior_specs]
    except ValueError as e:
        error_with_hint(str(e), ["Check --prior format: NAME=DIST:ARG1:ARG2"])
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
    try:
        likelihood_config = parse_likelihood(
            likelihood_spec,
            baseline_formula=baseline_formula,
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

    # --- Common settings ---
    method: str = getattr(args, "method", "mc")
    n_workers: int | None = getattr(args, "parallel", None)
    seed: int = getattr(args, "seed", 42)
    quiet: bool = getattr(args, "quiet", False)

    param_names = [p.name for p in priors]

    # --- Print summary ---
    if not quiet:
        print(f"=== tidal sample ({method}) ===")
        print(f"  Spec: {spec_path}")
        print(f"  Parameters: {', '.join(param_names)}")
        for p in priors:
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
        n_samples = getattr(args, "n_samples", 100)
        if not quiet:
            print(f"  Samples: {n_samples}")
            if n_workers:
                print(f"  Workers: {n_workers}")
            print()

        from tidal.inference._mc import run_monte_carlo

        result = run_monte_carlo(
            priors=priors,
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

            nlive = recommend_nlive(len(priors), nlive_auto)
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

        result = run_nested_sampling(
            log_likelihood=likelihood_fn,
            prior_transform=build_prior_transform(priors),
            ndim=len(priors),
            param_names=param_names,
            sampler="polychord",
            nlive=nlive,
            n_workers=n_workers,
            quiet=quiet,
            **ns_kwargs,
        )
        # Record priors in metadata so post-hoc analysis (e.g. correct
        # marginal D_KL per parameter, #308) can transform each column
        # into the space where its prior is uniform.
        result.metadata["priors"] = [
            {
                "name": p.name,
                "distribution": p.distribution,
                "low": p.low,
                "high": p.high,
            }
            for p in priors
        ]

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
    # On HPC login nodes mpi4py imports fine but COMM_WORLD initialisation
    # aborts (PMI2 not available outside a SLURM allocation).  Catch any
    # error and treat as rank 0 — the only consequence is that all ranks
    # in a real MPI run would write the file (rank-0 logic is a wallclock
    # optimisation, not correctness-critical).
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
        from tidal.inference._visualize import plot_importance

        importance_path = output_path / "importance.png"
        try:
            imp = result.parameter_importance()
            plot_importance(imp, importance_path)
            if not quiet:
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
