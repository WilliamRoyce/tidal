"""
Example showing how to conditionally use TeX-based plotting.

This script demonstrates how to gracefully handle the case where TeX
is not installed while still providing high-quality plots when available.
"""

import matplotlib.pyplot as plt
import numpy as np

from torsion_gertsenshtein import has_tex_support
from torsion_gertsenshtein.plot_pgf import enable_pgf


def main() -> None:
    """Demonstrate conditional TeX usage."""
    # Check if TeX support is available
    if has_tex_support():
        print("TeX support detected - using high-quality LaTeX rendering")
        enable_pgf("xelatex")  # or "pdflatex"/"lualatex"
    else:
        print("TeX not available - using standard matplotlib fonts")
        # The enable_pgf function will automatically fall back
        enable_pgf("xelatex")  # This will warn and use fallback

    # Create a simple plot
    x = np.linspace(0, 2 * np.pi, 100)
    y = np.sin(x)

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label=r"$\sin(x)$")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.title("Example Plot with Conditional TeX Support")
    plt.legend()
    plt.grid(visible=True, alpha=0.3)
    plt.tight_layout()

    # Save in multiple formats
    plt.savefig("example_plot.png", dpi=300, bbox_inches="tight")

    if has_tex_support():
        # Only save PDF if TeX is available for vector graphics
        plt.savefig("example_plot.pdf", bbox_inches="tight")
        print("Saved both PNG and PDF versions")
    else:
        print("Saved PNG version (PDF requires TeX)")

    plt.show()


if __name__ == "__main__":
    main()
