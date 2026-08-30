# Literature — Local arXiv TeX Sources

TeX source files for frequently referenced papers, kept locally so they can be read
directly instead of re-fetched from the web (see `.claude/rules/literature.md`: check here
before searching online, and read the TeX rather than the PDF).

**`docs/references.md` is the canonical index** — it is curated, grouped by topic, and
records *why* each paper is in the library. This file is the on-disk **inventory**: which
directories exist and which file to open in each.

The sources themselves are git-ignored (`.gitignore` → `literature/*`); this README is the
single tracked file in the directory, so the inventory travels with the repo while the
tarballs do not.

## Papers

The table is **generated from the directories on disk — do not hand-edit it.** It drifted
badly while hand-maintained (15 of 124 rows, and a header claiming 20). Regenerate after
adding or removing a paper:

```bash
uv run python -m scripts.bibaudit.index_literature          # rewrite the table
uv run python -m scripts.bibaudit.index_literature --check  # exit 1 if stale
```

Titles and authors come from the arXiv API (cached under `scripts/bibaudit/cache/`, so
re-runs need no network); the main TeX is the file containing `\begin{document}`. A few
pre-arXiv or publisher-PDF entries have no TeX and are marked `-- (PDF only)`.

<!-- BEGIN generated index -->

| Directory | Source | Main TeX | Title | Authors |
| --- | --- | --- | --- | --- |
| `0711.4866/` | [0711.4866](https://arxiv.org/abs/0711.4866) | `secDM5.tex` | Secluded WIMP Dark Matter | Pospelov et al. |
| `0803.0862/` | [0803.0862](https://arxiv.org/abs/0803.0862) | `xPerm.tex` | xPerm: fast index canonicalization for tensor computer algebra | Martin-Garcia |
| `0803.1967/` | [0803.1967](https://arxiv.org/abs/0803.1967) | `rmagnus.tex` | Describing neutrino oscillations in matter with Magnus expansion | Ioannisian & Smirnov |
| `0804.4011/` | [0804.4011](https://arxiv.org/abs/0804.4011) | `main.tex` | On Possible Light-Torsion Mixing in Background Magnetic Field | Kruglov |
| `0810.2989/` | [0810.2989](https://arxiv.org/abs/0810.2989) | `main.tex` | Confining potential from interacting magnetic and torsion fields | Gaete & Helaÿel-Neto |
| `0810.5488/` | [0810.5488](https://arxiv.org/abs/0810.5488) | `magnus-temp.tex` | The Magnus expansion and some of its applications | Blanes et al. |
| `0811.1030/` | [0811.1030](https://arxiv.org/abs/0811.1030) | `VbelowZ.tex` | Secluded U(1) below the weak scale | Pospelov |
| `0905.3732/` | [0905.3732](https://arxiv.org/abs/0905.3732) | `main.tex` | Infrared Modified Gravity with Dynamical Torsion | Nikiforova et al. |
| `0908.0629/` | [0908.0629](https://arxiv.org/abs/0908.0629) | `0908.0629.tex` | Constraints on background torsion from birefringence of CMB polarization | Das et al. |
| `0912.3767/` | [0912.3767](https://arxiv.org/abs/0912.3767) | `main.tex` | On the interplay between screening and confinement from interacting electromagnetic and torsion fields | Gaete & Helaÿel-Neto |
| `1104.2933/` | [1104.2933](https://arxiv.org/abs/1104.2933) | `approximations2.tex` | The Cosmic Linear Anisotropy Solving System (CLASS) II: Approximation schemes | Blas et al. |
| `1211.0500/` | [1211.0500](https://arxiv.org/abs/1211.0500) | `conversion-gw-gamma10.tex` | Conversion of relic gravitational waves into photons in cosmological magnetic fields | Dolgov & Ejlli |
| `1302.3884/` | [1302.3884](https://arxiv.org/abs/1302.3884) | `DP_stars_revision.tex` | New stellar constraints on dark photons | An et al. |
| `1308.3493/` | [1308.3493](https://arxiv.org/abs/1308.3493) | `xTras.tex` | xTras: a field-theory inspired xAct package for Mathematica | Nutma |
| `1401.4173/` | [1401.4173](https://arxiv.org/abs/1401.4173) | `article_2.tex` | Massive Gravity | Rham |
| `1405.7004/` | [1405.7004](https://arxiv.org/abs/1405.7004) | `cmb-lensing-prd.tex` | Effects of modified gravity on B-mode polarization | Amendola et al. |
| `1406.6646/` | [1406.6646](https://arxiv.org/abs/1406.6646) | `main.tex` | Canonical variational completion of differential equations | Voicu & Krupka |
| `1408.3978/` | [1408.3978](https://arxiv.org/abs/1408.3978) | `june2014revised.tex` | Advanced Virgo: a 2nd generation interferometric gravitational wave detector | Acernese et al. |
| `1411.4547/` | [1411.4547](https://arxiv.org/abs/1411.4547) | -- (PDF only) | Advanced LIGO | The LIGO Scientific Collaboration |
| `1506.02210/` | [1506.02210](https://arxiv.org/abs/1506.02210) | `main.tex` | The Theorem of Ostrogradsky | Woodard |
| `1508.02401/` | [1508.02401](https://arxiv.org/abs/1508.02401) | `curvaturesquaregravityarxiv.tex` | A Stueckelberg Approach to Quadratic Curvature Gravity and its Decoupling Limits | Hinterbichler & Saravani |
| `1510.06699/` | [1510.06699](https://arxiv.org/abs/1510.06699) | `ewgt-paper-submit2.tex` | Scale-invariant gauge theories of gravity: theoretical foundations | Lasenby & Hobson |
| `1702.05185/` | [1702.05185](https://arxiv.org/abs/1702.05185) | `main.tex` | Gravitational waves in Poincaré gauge gravity theory | Obukhov |
| `1706.07421/` | [1706.07421](https://arxiv.org/abs/1706.07421) | `higherderivpaper.tex` | Fixing extensions to General Relativity in the non-linear regime | Cayuso et al. |
| `1710.01562/` | [1710.01562](https://arxiv.org/abs/1710.01562) | `main.tex` | Perturbative reduction of derivative order in EFT | Glavan |
| `1804.05556/` | [1804.05556](https://arxiv.org/abs/1804.05556) | `main.tex` | General Poincaré gauge theory: Hamiltonian structure and particle spectrum | Blagojević & Cvetković |
| `1806.11020/` | [1806.11020](https://arxiv.org/abs/1806.11020) | `conversion_cqg_vf.tex` | Graviton-photon oscillation in alternative theories of gravity | Cembranos et al. |
| `1807.00171/` | [1807.00171](https://arxiv.org/abs/1807.00171) | `GRAPH_mixing.tex` | GRAPH mixing | Ejlli & Thandlam |
| `1807.10168/` | [1807.10168](https://arxiv.org/abs/1807.10168) | `template.tex` | Projective symmetries and induced electromagnetism in metric-affine gravity | Janssen & Jimenez-Cano |
| `1809.11095/` | [1809.11095](https://arxiv.org/abs/1809.11095) | `Paper.tex` | Rapid numerical solutions for the Mukhanov-Sasaki equation | Haddadin & Handley |
| `1812.02675/` | [1812.02675](https://arxiv.org/abs/1812.02675) | `PGT-particles.tex` | Ghost and tachyon free Poincaré gauge theories: a systematic approach | Lin et al. |
| `1903.02263/` | [1903.02263](https://arxiv.org/abs/1903.02263) | `main.tex` | Entropy in Poincaré gauge theory: Hamiltonian approach | Blagojević & Cvetković |
| `1905.04768/` | [1905.04768](https://arxiv.org/abs/1905.04768) | `paper.tex` | anesthetic: nested sampling visualisation | Handley |
| `1906.01421/` | [1906.01421](https://arxiv.org/abs/1906.01421) | `main.tex` | An efficient method for solving highly oscillatory ordinary differential equations with applications to physical systems | Agocs et al. |
| `1908.00232/` | [1908.00232](https://arxiv.org/abs/1908.00232) | `UHF_GW_Upper_limit_v3.tex` | Upper limits on the amplitude of ultra-high-frequency gravitational waves from graviton-photon mixing | Ejlli et al. |
| `1910.07506/` | [1910.07506](https://arxiv.org/abs/1910.07506) | `StabilityPGTsEJPC.tex` | Revisiting the Stability of Quadratic Poincaré Gauge Gravity | Jiménez & Torralba |
| `2003.00664/` | [2003.00664](https://arxiv.org/abs/2003.00664) | `Weyl_inv_cosmo0427.tex` | Consistent inflationary cosmology from quadratic gravity with dynamical torsion | Aoki & Mukohyama |
| `2003.02690/` | [2003.02690](https://arxiv.org/abs/2003.02690) | `prd_paper.tex` | Addressing $H_0$ tension with emergent dark radiation in unitary gravity | Barker et al. |
| `2004.02714/` | [2004.02714](https://arxiv.org/abs/2004.02714) | `main.tex` | Graviton-photon mixing. Exact solution in a constant magnetic field | Ejlli |
| `2005.01515/` | [2005.01515](https://arxiv.org/abs/2005.01515) | `master.tex` | The Dark Photon | Fabbrichesi et al. |
| `2005.02228/` | [2005.02228](https://arxiv.org/abs/2005.02228) | `WGT-particles.tex` | Ghost and tachyon free Weyl gauge theories: a systematic approach | Lin et al. |
| `2005.05574/` | [2005.05574](https://arxiv.org/abs/2005.05574) | `ptep01_main.tex` | Overview of KAGRA: Detector design and construction history | Akutsu et al. |
| `2006.01161/` | [2006.01161](https://arxiv.org/abs/2006.01161) | `text.tex` | Potential of radio telescopes as high-frequency gravitational wave detectors | Domcke & Garcia-Cely |
| `2006.03581/` | [2006.03581](https://arxiv.org/abs/2006.03581) | `apstemplate.tex` | Mapping Poincaré gauge cosmology to Horndeski theory for emergent dark energy | Barker et al. |
| `2008.09053/` | [2008.09053](https://arxiv.org/abs/2008.09053) | `source.tex` | Fresh perspective on gauging the conformal group | Hobson & Lasenby |
| `2009.05459/` | [2009.05459](https://arxiv.org/abs/2009.05459) | `varcompl_EPJP_format_Rev1_arxFin.tex` | Canonical variational completion and 4D Gauss-Bonnet gravity | Hohmann et al. |
| `2009.11739/` | [2009.11739](https://arxiv.org/abs/2009.11739) | `PGTandMG1120.tex` | Non-linearly ghost-free higher curvature gravity | Aoki |
| `2011.11254/` | [2011.11254](https://arxiv.org/abs/2011.11254) | `main.tex` | New Extraction of the Cosmic Birefringence from the Planck 2018 Polarization Data | Minami & Komatsu |
| `2011.12414/` | [2011.12414](https://arxiv.org/abs/2011.12414) | `Manuscript.tex` | Challenges and Opportunities of Gravitational Wave Searches at MHz to GHz Frequencies | Aggarwal et al. |
| `2101.02645/` | [2101.02645](https://arxiv.org/abs/2101.02645) | `apstemplate.tex` | Nonlinear Hamiltonian analysis of new quadratic torsion theories Part I. Cases with curvature-free constraints | Barker et al. |
| `2102.05859/` | [2102.05859](https://arxiv.org/abs/2102.05859) | `detection_paper.tex` | Rare Events Detected with a Bulk Acoustic Wave High Frequency Gravitational Wave Antenna | Goryachev et al. |
| `2102.10579/` | [2102.10579](https://arxiv.org/abs/2102.10579) | `main.tex` | General method for including Stueckelberg fields | Lyakhovich |
| `2105.04565/` | [2105.04565](https://arxiv.org/abs/2105.04565) | `main.tex` | Dark photon limits: a handbook | Caputo et al. |
| `2106.09355/` | [2106.09355](https://arxiv.org/abs/2106.09355) | `main.tex` | Reducible Stueckelberg symmetry and dualities | Abakumova & Lyakhovich |
| `2109.09546/` | [2109.09546](https://arxiv.org/abs/2109.09546) | `main.tex` | Gravity with dynamical torsion | Katanaev |
| `2110.13319/` | [2110.13319](https://arxiv.org/abs/2110.13319) | `GoryachevComment.tex` | Did Goryachev et al. detect megahertz gravitational waves? | Lasky & Thrane |
| `2111.14199/` | [2111.14199](https://arxiv.org/abs/2111.14199) | `draft_rot_ps.tex` | Computing Microwave Background Polarization Power Spectra from Cosmic Birefringence | Cai & Guan |
| `2112.11465/` | [2112.11465](https://arxiv.org/abs/2112.11465) | `main_v2.tex` | Detecting High-Frequency Gravitational Waves with Microwave Cavities | Berlin et al. |
| `2202.00032/` | [2202.00032](https://arxiv.org/abs/2202.00032) | `Ver2_MNRAS_Oct2023.tex` | Gertsenshtein-Zel$'$dovich effect: A plausible explanation for fast radio bursts? | Kushwaha et al. |
| `2202.00695/` | [2202.00695](https://arxiv.org/abs/2202.00695) | `GWaxion.tex` | A novel search for high-frequency gravitational waves with low-mass axion haloscopes | Domcke et al. |
| `2202.13919/` | [2202.13919](https://arxiv.org/abs/2202.13919) | `main_v5.tex` | New physics from the polarised light of the cosmic microwave background | Komatsu |
| `2204.06302/` | [2204.06302](https://arxiv.org/abs/2204.06302) | `Paper_Heating2018.tex` | Constraints on Primordial Magnetic Fields from their impact on the ionization history with Planck 2018 | Paoletti et al. |
| `2205.13534/` | [2205.13534](https://arxiv.org/abs/2205.13534) | `apstemplate.tex` | Geometric multipliers and partial teleparallelism in Poincaré gauge theory | Barker |
| `2205.13962/` | [2205.13962](https://arxiv.org/abs/2205.13962) | `ms.tex` | Improved Constraints on Cosmic Birefringence from the WMAP and Planck Cosmic Microwave Background Polarization Data | Eskilt & Komatsu |
| `2206.00658/` | [2206.00658](https://arxiv.org/abs/2206.00658) | `apstemplate.tex` | Supercomputers against strong coupling in gravity with curvature and torsion | Barker |
| `2208.03011/` | [2208.03011](https://arxiv.org/abs/2208.03011) | `main.tex` | Comparing Equivalent Gravities: common features and differences | Capozziello et al. |
| `2209.07804/` | [2209.07804](https://arxiv.org/abs/2209.07804) | `main.tex` | Isotropic cosmic birefringence from early dark energy | Murai et al. |
| `2210.15980/` | [2210.15980](https://arxiv.org/abs/2210.15980) | `main.tex` | A discrete discontinuity between the two phases of gravity | Sengupta |
| `2212.06924/` | [2212.06924](https://arxiv.org/abs/2212.06924) | `main-expanded.tex` | An adaptive spectral method for oscillatory second-order linear ODEs with frequency-independent cost | Agocs & Barnett |
| `2301.02072/` | [2301.02072](https://arxiv.org/abs/2301.02072) | `main.tex` | A Simple Derivation of the Gertsenshtein Effect | Palessandro & Rothman |
| `2302.03545/` | [2302.03545](https://arxiv.org/abs/2302.03545) | `main.tex` | The effective field theory approach to the strong coupling issue in $f(T)$ gravity | Hu et al. |
| `2302.08186/` | [2302.08186](https://arxiv.org/abs/2302.08186) | `conversionDef4arxiv.tex` | Graviton-photon oscillation in a cosmic background for a general theory of gravity | Cembranos et al. |
| `2303.11094/` | [2303.11094](https://arxiv.org/abs/2303.11094) | `Manuscript.tex` | Does gravitational confinement sustain flat galactic rotation curves without dark matter? | Barker et al. |
| `2308.09178/` | [2308.09178](https://arxiv.org/abs/2308.09178) | `article.tex` | Testing gravity with gauge-invariant polarization states of gravitational waves: Theory and pulsar timing sensitivity | Alves |
| `2309.14783/` | [2309.14783](https://arxiv.org/abs/2309.14783) | `source.tex` | Manifestly covariant variational principle for gauge theories of gravity | Hobson et al. |
| `2310.04150/` | [2310.04150](https://arxiv.org/abs/2310.04150) | `main.tex` | On graviton-photon conversions in magnetic environments | Hwang & Noh |
| `2311.03291/` | [2311.03291](https://arxiv.org/abs/2311.03291) | `main.tex` | DISCO-DJ I: a differentiable Einstein-Boltzmann solver for cosmology | Hahn et al. |
| `2311.11790/` | [2311.11790](https://arxiv.org/abs/2311.11790) | `Manuscript.tex` | Particle spectra of gravity based on internal symmetry of quantum fields | Barker |
| `2311.17459/` | [2311.17459](https://arxiv.org/abs/2311.17459) | `Hamiltonian_metric-affine_R2.tex` | Hamiltonian analysis of metric-affine-$R^2$ theory | Glavan et al. |
| `2312.17636/` | [2312.17636](https://arxiv.org/abs/2312.17636) | `main.tex` | Inverse Gertsenshtein effect as a probe of high-frequency gravitational waves | He et al. |
| `2401.15965/` | [2401.15965](https://arxiv.org/abs/2401.15965) | `GWPH-Longv8-Arxiv.tex` | Resonant Graviton-Photon Conversion with Stochastic Magnetic Field in the Expanding Universe | Addazi et al. |
| `2402.07641/` | [2402.07641](https://arxiv.org/abs/2402.07641) | `main.tex` | Particle spectra of general Ricci-type Palatini or metric-affine theories | Barker & Marzo |
| `2402.14917/` | [2402.14917](https://arxiv.org/abs/2402.14917) | `Manuscript.tex` | Consistent particle physics in metric-affine gravity from extended projective symmetry | Barker & Zell |
| `2403.15564/` | [2403.15564](https://arxiv.org/abs/2403.15564) | `main.tex` | Metric-affine cosmological models and the inverse problem of the calculus of variations. Part 1: variational bootstrapping -- the method | Ducobu & Voicu |
| `2405.01407/` | [2405.01407](https://arxiv.org/abs/2405.01407) | `GP8.tex` | Graviton-Photon Oscillations as a Probe of Quantum Gravity | Palessandro |
| `2405.08865/` | [2405.08865](https://arxiv.org/abs/2405.08865) | `main.tex` | Numerical Analysis of Resonant Axion-Photon Mixing: Part I | Ginés et al. |
| `2405.11786/` | [2405.11786](https://arxiv.org/abs/2405.11786) | `main.tex` | Graviton-photon conversions in Euler-Heisenberg nonlinear electrodynamics | Hwang & Noh |
| `2405.15581/` | [2405.15581](https://arxiv.org/abs/2405.15581) | `main.tex` | Nonlinear studies of modifications to general relativity: Comparing different approaches | Corman et al. |
| `2406.09500/` | [2406.09500](https://arxiv.org/abs/2406.09500) | `Manuscript.tex` | PSALTer: Particle Spectrum for Any Tensor Lagrangian | Barker et al. |
| `2406.09540/` | [2406.09540](https://arxiv.org/abs/2406.09540) | `part2_R1.tex` | Metric-affine cosmological models and the inverse problem of the calculus of variations. Part II: Variational bootstrapping of the $Λ$CDM model | Ducobu & Voicu |
| `2406.11956/` | [2406.11956](https://arxiv.org/abs/2406.11956) | `gravity_strong-CP_final_subm.tex` | Weyl-invariant Einstein-Cartan gravity: unifying the strong CP and hierarchy puzzles | Karananas et al. |
| `2406.12826/` | [2406.12826](https://arxiv.org/abs/2406.12826) | `Manuscript.tex` | Every Poincaré gauge theory is conformal: a compelling case for dynamical vector torsion | Barker et al. |
| `2406.17853/` | [2406.17853](https://arxiv.org/abs/2406.17853) | `main.tex` | Constraining gravitational-wave backgrounds from conversions into photons in the Galactic magnetic field | Lella et al. |
| `2406.18634/` | [2406.18634](https://arxiv.org/abs/2406.18634) | `main.tex` | Resonant Conversion of Gravitational Waves in Neutron Star Magnetospheres | McDonald & Ellis |
| `2407.09598/` | [2407.09598](https://arxiv.org/abs/2407.09598) | `curv_squared.tex` | The particle content of $R^2$ gravity revisited | Karananas |
| `2408.16818/` | [2408.16818](https://arxiv.org/abs/2408.16818) | `arbitrary_curv_squared.tex` | The particle content of (scalar curvature)$^2$ metric-affine gravity | Karananas |
| `2410.01355/` | [2410.01355](https://arxiv.org/abs/2410.01355) | `Torsion06.tex` | Classical electrodynamics and spin-torsion coupling effects | Trukhanova & Obukhov |
| `2410.03422/` | [2410.03422](https://arxiv.org/abs/2410.03422) | `main.tex` | Lowering the strong coupling mode of modified teleparallel gravity theories | Hu et al. |
| `2411.16928/` | [2411.16928](https://arxiv.org/abs/2411.16928) | `Stueckelberg_v2.tex` | Tensor global symmetries and the Stueckelberg mechanism for tensor fields | Chatzistavrakidis et al. |
| `2501.11723/` | [2501.11723](https://arxiv.org/abs/2501.11723) | `main.tex` | Challenges and Opportunities of Gravitational Wave Searches above 10 kHz | Aggarwal et al. |
| `2502.12517/` | [2502.12517](https://arxiv.org/abs/2502.12517) | `arXiv-Ver-3.tex` | Constraining circular polarization of high-frequency gravitational waves with CMB | Kushwaha & Jain |
| `2503.18972/` | [2503.18972](https://arxiv.org/abs/2503.18972) | `f_R,Q_.tex` | Vulnerability of f(Q) gravity theory and a possible resolution | Saha & Sanyal |
| `2505.00082/` | [2505.00082](https://arxiv.org/abs/2505.00082) | `Paper.tex` | Stable non-linear evolution in regularised higher derivative effective field theories | Figueras et al. |
| `2505.08457/` | [2505.08457](https://arxiv.org/abs/2505.08457) | `main.tex` | Graviton-photon conversion in blazar jets as a probe of high-frequency gravitational waves | Matsuo & Ito |
| `2505.23894/` | [2505.23894](https://arxiv.org/abs/2505.23894) | `Manuscript.tex` | Can metric-affine gravity be saved? | Barker et al. |
| `2506.02111/` | [2506.02111](https://arxiv.org/abs/2506.02111) | `Manuscript.tex` | The particle spectra of parity-violating theories: A less radical approach and an upgrade of PSALTer | Barker et al. |
| `2506.03609/` | [2506.03609](https://arxiv.org/abs/2506.03609) | `PBH.tex` | Experimental Limits on Planetary Mass Primordial Black Hole Mergers | Campbell et al. |
| `2506.07305/` | [2506.07305](https://arxiv.org/abs/2506.07305) | `r2v6.tex` | Hamiltonian equations of motion of quadratic gravity | Bellorin |
| `2506.17017/` | [2506.17017](https://arxiv.org/abs/2506.17017) | `DEMAG.tex` | Cosmology of Cubic Poincaré Gauge gravity | Bahamonde et al. |
| `2506.21662/` | [2506.21662](https://arxiv.org/abs/2506.21662) | `Manuscript.tex` | Infrared foundations for quantum geometry I: Catalogue of totally symmetric rank-three field theories | Barker et al. |
| `2507.02362/` | [2507.02362](https://arxiv.org/abs/2507.02362) | `main.tex` | Coupling Electromagnetism to Torsion: Black Holes and Spin-Charge Interactions | Bahamonde et al. |
| `2507.05349/` | [2507.05349](https://arxiv.org/abs/2507.05349) | `Manuscript.tex` | Infrared foundations for quantum geometry II: Catalogue of all torsion-like theories including new ghost-tachyon-free cases | Barker et al. |
| `2507.09228/` | [2507.09228](https://arxiv.org/abs/2507.09228) | `paper_Qtorsion.tex` | Alleviating the Hubble tension with Torsion Condensation (TorC) | Legner et al. |
| `2507.16609/` | [2507.16609](https://arxiv.org/abs/2507.16609) | `Gert.tex` | Gravitational Wave Scattering on Magnetic Fields | Domcke et al. |
| `2510.08201/` | [2510.08201](https://arxiv.org/abs/2510.08201) | `spectrum_R2.tex` | Spectrum of pure $R^2$ gravity: full Hamiltonian analysis | Barker & Glavan |
| `2510.17094/` | [2510.17094](https://arxiv.org/abs/2510.17094) | `Final_version.tex` | Gertsenshtein effect on the spacetime curved by background magnetic field with geometric optics | Tomomatsu et al. |
| `2512.25007/` | [2512.25007](https://arxiv.org/abs/2512.25007) | `Manuscript.tex` | Fast Poisson brackets and constraint algebras in canonical gravity | Barker |
| `2601.22007/` | [2601.22007](https://arxiv.org/abs/2601.22007) | `main.tex` | Stückelberg inspired approach for avoiding singular Hamiltonians in Lorentz violating models of antisymmetric tensor field | Aashish & Saif |
| `2602.12114/` | [2602.12114](https://arxiv.org/abs/2602.12114) | `draft_sfj.tex` | Matrix bordering structure of the Faddeev-Jackiw algorithm: kernel reduction and symbolic automation | Chan-López et al. |
| `2602.23466/` | [2602.23466](https://arxiv.org/abs/2602.23466) | `main.tex` | nanoCMB: A minimal CMB power spectrum calculator in Python | Moss |
| `2604.12775/` | [2604.12775](https://arxiv.org/abs/2604.12775) | `GZ20260413.tex` | Gravitational Gertsenshtein-Zeldovich mechanism for the Association between GW190425 and FRB 20190425A | Wu et al. |
| `2606.30785/` | [2606.30785](https://arxiv.org/abs/2606.30785) | `Manuscript.tex` | Numerical polology: towards next-generation model-building for cosmology | Barker et al. |
| `2608.06480/` | [2608.06480](https://arxiv.org/abs/2608.06480) | `main.tex` | Cosmic birefringence from a joint analysis of ACT and Planck | Eskilt |
| `BF02721794/` | [BF02721794](https://doi.org/10.1007/BF02721794) | -- (PDF only) | Hamiltonian structure of the theory of gravity with R+T^2 type of Lagrangian | Blagojevic & Nikolic (1983) |
| `PhysRep258.1/` | [PhysRep258.1](https://doi.org/10.1016/0370-1573(94)00111-F) | -- (PDF only) | Metric-affine gauge theory of gravity: field equations, Noether identities, world spinors, and breaking of dilation invariance | Hehl, McCrea, Mielke & Ne'eman (1995) |
| `PhysRevD.28.2455/` | [PhysRevD.28.2455](https://doi.org/10.1103/PhysRevD.28.2455) | -- (PDF only) | Hamiltonian dynamics of Poincare gauge theory: General structure in the time gauge | Blagojevic & Nikolic (1983) |
| `astro-ph_9506072/` | [astro-ph/9506072](https://arxiv.org/abs/astro-ph/9506072) | `9506072.tex` | Cosmological Perturbation Theory in the Synchronous and Conformal Newtonian Gauges | Ma & Bertschinger |
| `astro-ph_9603033/` | [astro-ph/9603033](https://arxiv.org/abs/astro-ph/9603033) | `los.tex` | A Line of Sight Approach to Cosmic Microwave Background Anisotropies | Seljak & Zaldarriaga |
| `gr-qc_0001010/` | [gr-qc/0001010](https://arxiv.org/abs/gr-qc/0001010) | `gyros17.tex` | How does the electromagnetic field couple to gravity, in particular to metric, nonmetricity, torsion, and curvature? | Hehl & Obukhov |
| `gr-qc_0112030/` | [gr-qc/0112030](https://arxiv.org/abs/gr-qc/0112030) | `main.tex` | Hamiltonian Analysis of Poincaré Gauge Theory: Higher Spin Modes | Yo & Nester |
| `gr-qc_0302040/` | [gr-qc/0302040](https://arxiv.org/abs/gr-qc/0302040) | `kop1.tex` | Three lectures on Poincare gauge theory | Blagojevic |
| `gr-qc_0305049/` | [gr-qc/0305049](https://arxiv.org/abs/gr-qc/0305049) | `main.tex` | Torsion nonminimally coupled to the electromagnetic field and birefringence | Rubilar et al. |
| `gr-qc_0307063/` | [gr-qc/0307063](https://arxiv.org/abs/gr-qc/0307063) | `main.tex` | Maxwell's field coupled nonminimally to quadratic torsion: Induced axion field and birefringence of the vacuum | Itin & Hehl |
| `gr-qc_0311082/` | [gr-qc/0311082](https://arxiv.org/abs/gr-qc/0311082) | `GRET-jhep.tex` | Quantum Gravity in Everyday Life: General Relativity as an Effective Field Theory | Burgess |
| `gr-qc_9211002/` | [gr-qc/9211002](https://arxiv.org/abs/gr-qc/9211002) | `main.tex` | Einstein Equation with Quantum Corrections Reduced to Second Order | Parker & Simon |
| `gr-qc_9402012/` | [gr-qc/9402012](https://arxiv.org/abs/gr-qc/9402012) | `main.tex` | Metric-Affine Gauge Theory of Gravity: Field Equations, Noether Identities, World Spinors, and Breaking of Dilation Invariance | Hehl et al. |
| `gr-qc_9405057/` | [gr-qc/9405057](https://arxiv.org/abs/gr-qc/9405057) | `9405057.tex` | General relativity as an effective field theory: The leading quantum corrections | Donoghue |
| `gr-qc_9902032/` | [gr-qc/9902032](https://arxiv.org/abs/gr-qc/9902032) | `main.tex` | Hamiltonian analysis of Poincaré gauge theory scalar modes | Yo & Nester |
| `hep-ph_9306321/` | [hep-ph/9306321](https://arxiv.org/abs/hep-ph/9306321) | `main.tex` | Effective Lagrangians with Higher Order Derivatives | Grosse-Knetter |
| `hep-th_0103093/` | [hep-th/0103093](https://arxiv.org/abs/hep-th/0103093) | `main.tex` | Physical Aspects of the Space-Time Torsion | Shapiro |
| `hep-th_0301256/` | [hep-th/0301256](https://arxiv.org/abs/hep-th/0301256) | `main.tex` | Irregular Hamiltonian Systems | Miskovic & Zanelli |
| `hep-th_0302033/` | [hep-th/0302033](https://arxiv.org/abs/hep-th/0302033) | `main.tex` | Dynamical Structure of Irregular Constrained Systems | Miskovic & Zanelli |
| `hep-th_0406216/` | [hep-th/0406216](https://arxiv.org/abs/hep-th/0406216) | `koganweb.tex` | Heisenberg-Euler Effective Lagrangians : Basics and Extensions | Dunne |

<!-- END generated index -->

## Adding new papers

```bash
uv run python -m scripts.bibaudit.fetch_fulltext <arxiv-id> [<arxiv-id> ...]
```

This downloads each e-print, extracts the `.tex`/`.bbl` into `literature/<id>/` (old-style
ids such as `gr-qc/0305049` become `gr-qc_0305049`), skips anything already present, and
throttles between fetches. Then:

1. regenerate the table above, and
2. add a curated line to `docs/references.md` saying why the paper is in the library.
