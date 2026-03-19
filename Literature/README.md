# Literature — Local arXiv TeX Sources

TeX source files for frequently referenced papers, kept locally to avoid repeated web fetches.
This directory is git-ignored.

## Papers

| Directory | arXiv ID | Main TeX | Title | Authors |
|-----------|----------|----------|-------|---------|
| `2406.12826/` | [2406.12826](https://arxiv.org/abs/2406.12826) | `Manuscript.tex` | Every Poincare gauge theory is conformal: a compelling case for dynamical vector torsion | Barker et al. |
| `2301.02072/` | [2301.02072](https://arxiv.org/abs/2301.02072) | `main.tex` | A Simple Derivation of the Gertsenshtein Effect | Domcke & Garcia-Cely |
| `2310.04150/` | [2310.04150](https://arxiv.org/abs/2310.04150) | `main.tex` | On graviton-photon conversions in magnetic environments | Domcke & Garcia-Cely |
| `2406.17853/` | [2406.17853](https://arxiv.org/abs/2406.17853) | `main.tex` | Constraining GW backgrounds from conversions into photons in the Galactic magnetic field | Dandoy, Lella et al. |
| `2405.08865/` | [2405.08865](https://arxiv.org/abs/2405.08865) | `main.tex` | Numerical Analysis of Resonant Axion-Photon Mixing: Part I | Berlin et al. |
| `2507.16609/` | [2507.16609](https://arxiv.org/abs/2507.16609) | `Gert.tex` | Gravitational Wave Scattering on Magnetic Fields | Domcke, Garcia-Cely & Lee |
| `2410.01355/` | [2410.01355](https://arxiv.org/abs/2410.01355) | `Torsion06.tex` | Classical electrodynamics and spin-torsion coupling effects | Obukhov et al. |
| `hep-th_0103093/` | [hep-th/0103093](https://arxiv.org/abs/hep-th/0103093) | `main.tex` | Physical Aspects of the Space-Time Torsion | Shapiro |

## Adding new papers

Download TeX source from arXiv and extract:
```bash
mkdir -p Literature/<ID> && cd Literature/<ID>
curl -sL -o source.tar.gz "https://arxiv.org/e-print/<ID>"
tar xzf source.tar.gz 2>/dev/null || gunzip -c source.tar.gz > main.tex
rm source.tar.gz
```
Then update this README table.
