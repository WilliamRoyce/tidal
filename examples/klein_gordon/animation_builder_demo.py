"""Demonstration of AnimationBuilder for simplified animation creation.

This example shows how to use the new AnimationBuilder class to create
various types of animations with minimal boilerplate code.

The AnimationBuilder eliminates ~100-200 lines of repeated animation setup
code per example by providing a fluent, configuration-based API.
"""

from __future__ import annotations

from torsion_gertsenshtein.kgsim import (
    AnimationBuilder,
    AnimationConfig,
    GaussianPulse,
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
    run_with_snapshots,
)


def demo_1d_spacetime() -> None:
    """Demonstrate 1D spacetime plot creation."""
    print("\n" + "=" * 70)
    print("DEMO 1: 1D Spacetime Plot")
    print("=" * 70)

    # Setup and run simulation
    grid = make_grid(GridConfig(dim=1, shape=(128,), bounds=((0.0, 100.0),)))
    ic = GaussianPulse(amplitude=1.0, width=3.0, center=[50.0])
    state = ic.build(grid)

    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=50.0, solver="scipy", backend="numpy", progress=True
    )

    _result, storage = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=0.5
    )

    # Create animation using builder - just 3 lines!
    builder = AnimationBuilder(storage, grid)
    anim_config = AnimationConfig(
        output_path="outputs/spacetime_demo.png",
        title="1D Klein-Gordon Spacetime",
        xlabel="Position x",
        ylabel="Time t",
        cbar_label="φ(x,t)",
    )
    builder.create_spacetime_1d(anim_config)

    print("✓ Created spacetime plot with AnimationBuilder")
    print(f"  Snapshots: {len(storage)}")
    print("  Lines of animation code: ~10 (vs ~100 without builder)")


def demo_1d_line_animation() -> None:
    """Demonstrate 1D line animation creation."""
    print("\n" + "=" * 70)
    print("DEMO 2: 1D Line Animation")
    print("=" * 70)

    # Setup
    grid = make_grid(GridConfig(dim=1, shape=(128,), bounds=((0.0, 100.0),)))
    ic = GaussianPulse(amplitude=1.0, width=3.0, center=[30.0], initial_velocity=0.0)
    state = ic.build(grid)

    pde = KleinGordonPDE(params=KGParameters(mass=0.3))
    config = SimulationConfig(
        t_end=30.0, solver="scipy", backend="numpy", progress=True
    )

    _result, storage = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=0.3
    )

    # Create line animation
    builder = AnimationBuilder(storage, grid)
    anim_config = AnimationConfig(
        output_path="outputs/line_animation_demo",  # Extension auto-added
        title="1D Klein-Gordon Wave",
        xlabel="Position x",
        cbar_label="φ(x)",
        fps=20,
    )
    builder.create_1d_line_animation(anim_config)

    print("✓ Created 1D line animation")
    print(f"  Frames: {len(storage)}")


def demo_2d_heatmap_animation() -> None:
    """Demonstrate 2D heatmap animation creation."""
    print("\n" + "=" * 70)
    print("DEMO 3: 2D Heatmap Animation")
    print("=" * 70)

    # Setup
    grid = make_grid(
        GridConfig(dim=2, shape=(64, 64), bounds=((-20.0, 20.0), (-20.0, 20.0)))
    )
    ic = GaussianPulse(amplitude=1.0, width=3.0, center=[0.0, 0.0])
    state = ic.build(grid)

    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=20.0, solver="scipy", backend="numpy", progress=True
    )

    _result, storage = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=0.5
    )

    # Create 2D animation
    builder = AnimationBuilder(storage, grid)
    anim_config = AnimationConfig(
        output_path="outputs/2d_heatmap_demo",
        title="2D Klein-Gordon Evolution",
        xlabel="x",
        cmap="RdBu_r",
        figsize=(8, 8),
        fps=15,
    )
    builder.create_2d_heatmap_animation(anim_config)

    print("✓ Created 2D heatmap animation")
    print(f"  Frames: {len(storage)}")


def main() -> None:
    """Run all AnimationBuilder demonstrations."""
    print("\n" + "=" * 70)
    print("ANIMATION BUILDER DEMONSTRATION")
    print("=" * 70)
    print("""
The AnimationBuilder class provides:

✓ Unified API for all animation types
✓ Automatic writer selection (ffmpeg/pillow)
✓ Consistent colormap normalization
✓ Fluent configuration-based interface
✓ ~80% reduction in animation boilerplate

Old pattern (100-200 lines per animation):
    - Manual figure setup
    - Manual colormap normalization
    - Manual writer selection
    - Manual frame update functions
    - Duplicate code across examples

New pattern (5-15 lines per animation):
    - AnimationBuilder(storage, grid)
    - AnimationConfig(output_path, title, ...)
    - builder.create_XXX_animation(config)
    """)

    demo_1d_spacetime()
    demo_1d_line_animation()
    demo_2d_heatmap_animation()

    print("\n" + "=" * 70)
    print("CODE COMPARISON")
    print("=" * 70)
    print("""
BEFORE (old pattern for 2D animation):
    ```python
    # ~150 lines of code
    fig, ax = plt.subplots(...)
    # Manual colormap setup
    vmin, vmax = data.min(), data.max()
    if vmin < 0 < vmax:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    # Manual image setup
    im = ax.imshow(...)
    # Manual update function
    def update(frame):
        ...
    # Manual writer selection
    if shutil.which("ffmpeg"):
        writer = FFMpegWriter(...)
    else:
        writer = PillowWriter(...)
    # Create and save animation
    anim = FuncAnimation(...)
    anim.save(...)
    ```

AFTER (new AnimationBuilder pattern):
    ```python
    # ~10 lines of code
    builder = AnimationBuilder(storage, grid)
    config = AnimationConfig(
        output_path="output.mp4",
        title="Evolution",
        cmap="RdBu_r"
    )
    builder.create_2d_heatmap_animation(config)
    ```

Result: 93% less code, consistent API, easier to maintain!
    """)


if __name__ == "__main__":
    main()
