import matplotlib.pyplot as plt

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
)


def main() -> None:
    grid_config = GridConfig(
        dim=1, shape=(1024,), bounds=((0.0, 200.0),), periodic=True
    )
    grid = make_grid(grid_config)

    params = KGParameters(mass=0.5)
    pde = KleinGordonPDE(params)

    state = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,  # adaptive
        solver="scipy",  # or "explicit"
        method="RK45",
        backend="numba",
        progress=True,
    )

    result = run(pde=pde, state=state, config=simulation_config)

    # Plot φ at the end
    phi = result[0]
    phi.plot(title=r"Klein-Gordon: $\phi(x, t_{\mathrm{end}})$")

    plt.gcf().savefig("phi_end.png", dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
