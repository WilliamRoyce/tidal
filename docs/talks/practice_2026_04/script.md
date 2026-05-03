# Practice talk — spoken script

Continuous-prose script for live delivery. Written to be learned and rehearsed.
Total target ~14–16 min spoken, leaving ~5 min buffer in a 20-min slot.

Conventions:

- `▶ NEXT — …` marks each slide-advance click.
- Greek letters are written out phonetically (`kappa B-zero`, `alpha-i I-i`, `delta-one`) so the speaker doesn't trip on symbols.
- Prose is what comes out of the mouth; what's on the slide is referred to ("the formula on the slide", "the table") rather than read out.
- For each slide: target time + speaking notes + the `▶ NEXT` cue.

---

═══════════════════════════════════════════════════════════════
▶ Title slide (~15 s)
═══════════════════════════════════════════════════════════════

Hi everyone. I'm William, a master's student at the
Cavendish, supervised by Will Barker, Mike Hobson, and Anthony Lasenby, working on whether extensions to our theory of gravity can amplify a particular
old graviton-photon conversion process. This is a practice talk on my work the past few months on this project —
feedback on both content and delivery is very welcome.

▶ NEXT — slide 1: motivation

═══════════════════════════════════════════════════════════════
▶ Slide 1 — 01 motivation (~75 s)
═══════════════════════════════════════════════════════════════

The starting point is the Gertsenshtein effect: a graviton propagating
through a background magnetic field can, perhaps surprisingly, convert into a
photon. However, it is actually a textbook prediction of Einstein–Maxwell
theory, and it has been around since 1962.

Perhaps unfortunately, this effect is astrophysically
useless. For a magnetar with the strongest known magnetic fields,
the conversion probability is around ten-to-the-minus-ten. And for any
laboratory setup it's even more hopeless. At these scales, there is no foreseeable detector that will see Gertsenshtein
in any astrophysical or laboratory setting, within standard-GR.

Our theoretical curiosity however prompts us the question: is there any extension of GR that can amplify this conversion to a detectable level?

The approach
in this project is unapologetically phenomenological — instead of
picking one specific extension of GR and defending it, I'm
systematically scanning the space of linearised curvature, torsion, and electromagnetic
modifications and asking which subsectors, if any, amplify the
conversion.

▶ NEXT — slide 2: the Gertsenshtein effect

═══════════════════════════════════════════════════════════════
▶ Slide 2 — 02 gertsenshtein (~60 s)
═══════════════════════════════════════════════════════════════

So what actually is the Gertsenshtein effect. The setup on the slide:
linearised GR plus Maxwell in a uniform background magnetic field
B-zero. The reason the graviton and the photon mix is mechanical — the
EM action contains two metric factors, one in each F-down-mu-nu and
F-up-mu-nu. When you perturb the metric around B-zero, those two factors
generate a cross-term coupling the metric perturbation to the photon
perturbation.

What you see on the slide is the conversion-probability formula —
sin-squared in kappa B-zero D over two. The structure is exactly two-state
quantum-mechanical mixing, the same machinery you see in axion-photon
conversion through a Primakoff process. The intuition: the background
field acts like an off-diagonal coupling between the graviton and the
photon dispersion relations, with a beat frequency set by kappa B-zero.

So now the question becomes: how big is that beat frequency in any
realistic setting?

▶ NEXT — slide 3: why it's so weak

═══════════════════════════════════════════════════════════════
▶ Slide 3 — 03 literature (~50 s)
═══════════════════════════════════════════════════════════════

The answer to that question is on the slide: the conversion probability
scales as B-zero times distance over the Planck mass, squared. M-Planck
is ten-to-the-nineteen GeV — that's the suppression. Even for a magnetar
running over its full ten-kilometre scale, B-zero D divided by M-Planck
is around ten-to-the-minus-five, giving a probability of ten-to-the-minus-ten.

Three literature points worth knowing. Gertsenshtein himself in 1962
gave the formula and immediately recognised it was hopelessly small.
Boccaletti and collaborators in 1970 did the full analytic solution for
a localised B-field region. Raffelt and Stodolsky in 1988 embedded the
problem in the axion-photon-mixing framework, which is the language most
people use today.

The point to make explicit: this formula is exact in linearised
Einstein-Maxwell. There's no approximation to soften, no refinement that
brings it up. Any amplification has to come from new gravitational
physics outside Einstein-Maxwell.

▶ NEXT — slide 4: PGT and torsion

═══════════════════════════════════════════════════════════════
▶ Slide 4 — 04 pgt_a (~60 s)
═══════════════════════════════════════════════════════════════

The natural place to look is Poincaré gauge theory. The framing the
slide leads with: spacetime has two independent local symmetries —
local Lorentz invariance and local translations — and the slide's
table makes the analogy explicit. In QED you gauge U(1), get a
connection A-mu, and a field strength F. In gravity you can gauge
both subgroups of the Poincaré group: gauging the Lorentz part gives
you a spin connection and curvature, gauging translations gives you
the tetrad and torsion.

GR keeps one of these: curvature. PGT keeps both. So the bottom line:
in GR torsion is set to zero by assumption; in PGT it's a new
gravitational degree of freedom on equal footing with the metric.

The figure on the right is the geometric intuition. Parallel-transport a
vector around a closed loop and you get two mismatches: a rotational
mismatch — that's curvature — and a translational one where the loop
literally fails to close, which is torsion. The third panel is
non-metricity, which neither GR nor PGT keeps.

▶ NEXT — slide 5: the parameter landscape

═══════════════════════════════════════════════════════════════
▶ Slide 5 — 05 pgt_b (~60 s)
═══════════════════════════════════════════════════════════════

Once you commit to PGT, the parameter space gets large. The schematic
Lagrangian on the left has three families: torsion-squared invariants
giving torsion masses; curvature-squared invariants giving graviton
self-coupling and torsion mixing; and non-minimal couplings like
delta-one R-tilde-anti-symmetric F-mu-nu, which directly shuttle energy
between the gravitational and electromagnetic sectors. The
non-minimal terms are the most interesting because they're the only
ones that can directly modify the Gertsenshtein vertex — they're absent
in GR and they're the principal target of the scan.

The table on the right enumerates the sectors. Adding up the entries
gives roughly 150 independent quadratic couplings — and that's before
you add cubic vertices on the B-zero background, which is where the
Gertsenshtein-type physics actually lives. Any of these could in
principle amplify; hand-deriving any single one safely is fragile, and
hand-deriving all of them is impossible.

So the program is to scan this whole landscape systematically.

▶ NEXT — slide 6: the question and the method

═══════════════════════════════════════════════════════════════
▶ Slide 6 — 06 question (~45 s)
═══════════════════════════════════════════════════════════════

This is the central question, framed exactly as on the slide:
is there any stable combination of quadratic curvature and torsion
modifications whose linearised PDEs, in a background magnetic field,
give an amplified Gertsenshtein conversion?

For every Lagrangian we test, three things have to happen. We derive
the linearised field equations symbolically — by hand they're
intractable for anything but the simplest case. We screen for unstable
regimes by computing eigenvalues across all spatial modes — tachyonic
points get rejected before they enter the inference. And we compute the
conversion probability and compare it to the GR baseline.

All three of those steps are what TIDAL automates, and that's the
next slide.

▶ NEXT — slide 7: TIDAL overview

═══════════════════════════════════════════════════════════════
▶ Slide 7 — 07 tidal_overview (~60 s)
═══════════════════════════════════════════════════════════════

TIDAL stands for Tensor Integration and Derivation for Any Lagrangian.
The pipeline goes left-to-right across the slide: a TOML config file
declares the Lagrangian and the field content; Mathematica with the
xAct and xTras packages does the symbolic work — Euler-Lagrange
variation, multi-field perturbation around a background, and the
covariant-derivative expansion that brings in the contortion contributions; that produces a JSON spec; and the Python
solver consumes the JSON, runs sweeps, runs inference, generates plots.

The core claim, which I want to land properly: every equation in the
numerical solver is traceable back to a single Lagrangian declaration.
No hand-derived PDEs. No manual index gymnastics. For a phenomenological
scan that varies the Lagrangian across dozens of operator combinations,
the underlying tensor calculus and constraint reductions are
genuinely intractable by hand — symbolic computation isn't
a convenience here, it's the only viable route.

▶ NEXT — slide 8: the modal solver

═══════════════════════════════════════════════════════════════
▶ Slide 8 — 08 tidal_numerical (~70 s)
═══════════════════════════════════════════════════════════════

This is the numerical engine. It's a spectral solver tailored for the
linearised gauge theories the survey produces. Two regimes — per-mode
when the background magnetic field is uniform, and a convolution-coupled
form when it's localised in space — and both evolve by direct matrix
exponential. Crucially, no eigendecomposition.

The reason it fits this survey class: every theory we test linearises
to periodic, time-independent PDEs around a flat background, which are
exactly the modal solver's eligibility criteria. The benefit is
machine-precision time evolution with no time-stepping error to obscure
the weak signals we're trying to measure.

The thing worth flagging is what we had to do to make this work for
gauge theories. The textbook approach — V times diag-e-lambda-t times
V-inverse — fails universally on PGT because the eigenvector matrix V
is catastrophically ill-conditioned: condition numbers in the
ten-to-the-fifteen range, where inversion just doesn't work
numerically. So we compute the matrix exponential directly via the
Padé approximant, never forming V-inverse. Second fix: Fourier-space
algebraic elimination of fields that have no time derivative — the
Gauss-law-type constraints that PGT generates in abundance. Both
formulas are in the backup if anyone wants to see them.

▶ NEXT — slide 9: the inference framework

═══════════════════════════════════════════════════════════════
▶ Slide 9 — 09 inference (~75 s)
═══════════════════════════════════════════════════════════════

The inference layer turns the survey into a Bayesian search. For each
parameter point theta we run a full PDE simulation, get a peak
conversion probability, and divide it by the GR baseline — sin-squared
kappa B-zero t over two. That ratio is what we call the amplification
factor A of theta.

The slide shows the small-coupling limit: when the conversion is well
inside the perturbative regime, the amplification factor reduces to a
ratio of conversion coefficients independent of B-zero. That's the
quantity we actually scan over — it's a property of the Lagrangian,
not of the experimental setup. The likelihood for nested sampling is
plus-or-minus log A, depending on whether we're searching for
amplification or suppression.

Two summary statistics come out the back. First, the Bayesian evidence
Z, which equals the prior-averaged amplification factor — log Z near
zero means the prior is null on average, log Z much greater than zero
means amplification persists across the prior. Second, the
Kullback-Leibler divergence between the posterior and the prior, which
tells us whether the enhancement is broadly distributed across the
parameter space, meaning generic, or concentrated in a narrow corner,
meaning fine-tuned.

So when I report results I'm reporting these summary numbers, not raw
P-max values.

▶ NEXT — slide 10: pre-campaign results

═══════════════════════════════════════════════════════════════
▶ Slide 10 — 11 results (~90 s)
═══════════════════════════════════════════════════════════════

Three pre-campaign models that determined the structure of the full
HPC scan.

The first is minimal PGT: Einstein–Cartan, plus the three torsion-mass
invariants alpha-i I-i, plus b-five times R-tilde-squared. The plane-wave
reduced equations split into two algebraically decoupled channels, and
the Gertsenshtein channel — the one that an incoming TT graviton can
actually excite — contains zero torsion fields in either equation.
That's a structural property of the Lagrangian, visible directly in the
derived equations. So A equals one exactly, non-perturbatively in
b-five — no simulation required. This null is the reason the campaign
targets non-minimal sectors specifically.

The second is the plasma baseline: pure Einstein–Maxwell plus an
effective photon mass m-A-squared. The photon mass detunes the
graviton-photon resonance, and P-max gets monotonically suppressed.
This sets the floor any BSM amplification mechanism has to overcome —
plasma alone, in pure GR, suppresses the Gertsenshtein conversion.

The third is the dark-photon analogy. The slide shows the Lagrangian:
torsion trace as a hidden U(1) Proca vector, mass coming from the I-3
invariant on the full torsion tensor, kinetic mixing
delta-m F-dot-F-of-T into electromagnetism. In vacuum, an algebraic
cancellation makes the photon and the dark photon exact eigenmodes —
no graviton-induced channel opens. We then add plasma to break the
eigenmode degeneracy, which is the standard route from dark-photon
phenomenology — that's what gives you Raffelt-Stodolsky-type resonant
conversion. The result: kinetic mixing has negligible effect on total
conversion. Only the photon mass matters, and only as Gertsenshtein
detuning. A is at most one everywhere.

▶ NEXT — slide 11: campaign progress

═══════════════════════════════════════════════════════════════
▶ Slide 11 — 12 outlook (~75 s)
═══════════════════════════════════════════════════════════════

This is the campaign progress map. The table walks through stages A
through D2.1, all closed. Dark-photon plasma — null on amplification,
informative on suppression. Einstein–Cartan — null. R-squared PGT, the
structural one I described — null. Ricci-EM, the
delta-one R-tilde-mu-nu F-mu-nu non-minimal coupling — strong
suppression, Bayes factor of order ten-to-the-seven. Bahamonde
Yang-Mills PGT in five dimensions — null both ways. Barker chi-axial
in six dimensions — null both ways, the chi posterior centred at zero.

The story emerging: amplification is elusive. Multiple natural-looking
PGT extensions are Gertsenshtein-neutral. The campaign so far has been
informative by elimination — it's mapped out where amplification
_can't_ live.

What's queued: the Shapiro PGT and complete-PGT stages. The parity-odd
Yang-Mills extension — the abstract restricted to parity-even but parity-odd
is a natural follow-up. Plasma-Proca extensions where the photon's
effective mass comes from realistic plasma physics, relevant for
magnetar and IGM observability. And eventually the complete even-PGT
and odd-PGT landscapes with sixty-plus couplings. Each new theory is a
config-file change, not a code rewrite.

▶ NEXT — slide 12: summary

═══════════════════════════════════════════════════════════════
▶ Slide 12 — 13 summary (~45 s)
═══════════════════════════════════════════════════════════════

Three things to take away.

First, standard Gertsenshtein is astrophysically useless — any
detectable signal needs new gravitational physics. Second, we're
systematically scanning the quadratic curvature and torsion landscape,
with the inference cleanly separating amplification, suppression, and
fine-tuning. Third, the first half of the scan is done — simple PGT
extensions don't amplify, the non-minimal Ricci-EM direction gives
strong suppression, and the higher-dimensional, parity-odd, and
plasma corners are next.

The line on the slide is the one I want to land carefully:
whether or not amplification exists in this landscape, we're building
the systematic map of where it could live.

Thanks — questions and delivery feedback are both welcome.

▶ END — Q&A

═══════════════════════════════════════════════════════════════
Backup — for Q&A only
═══════════════════════════════════════════════════════════════

These are short notes for the appendix slides, used only if a question
takes me there. No need to memorise verbatim.

**B0 — Magnetar context.** Visual reference for "magnetar field
ten-to-the-fifteen Gauss over ten kilometres gives P around
ten-to-the-minus-ten". The slide is a single picture; just point at the
field strength and the scale.

**B1 — Higher-derivative terms via perturbative reduction.** This is
the slide to reach for if anyone asks how we actually included
R-tilde-squared in the minimal-PGT calculation, or about the
Ostrogradsky-ghost issue. The story: R-tilde-squared gives 4th-order
time derivatives, which generically have ghost modes. We use
Parker–Simon iterative reduction — treat b-five as small, expand
order by order, the 4th-order operator never acts on the unknown but
only as a known source built from the previous order's trajectory, so
the ghost branch never enters the propagator. The order-1 correction
is solved analytically by Duhamel's principle: convolution against the
base propagator, evaluated in closed form using the order-zero
eigenbasis. The constraint-promotion barrier — what stops us doing
energy measurements: at b-five equals zero, some metric components
are pure algebraic constraints; the R-tilde-squared correction
promotes them to 4th-order propagating fields, the number of
propagating modes jumps discontinuously, and there's no published
Hamiltonian-side recipe that handles this case. That's why the
minimal-PGT null result rests on a structural argument — the
polarisation block-diagonalisation visible directly in the EOMs —
rather than on an energy comparison.

**B2 — Modal solver: two structural fixes.** The slide to reach for
if anyone asks about the solver. Per-mode k satisfies y-dot equals
M-of-k times y; we evaluate y-of-t as exp-M-of-k-times-t times y-zero.
Two fixes: (i) ill-conditioned eigenvectors — cond V is ten-to-the-fifteen
to ten-to-the-seventeen for gauge theories, so we use Padé directly
instead of V-diag-V-inverse; (ii) constraint fields, which enter
algebraically — partition into dynamic and constrained, eliminate the
constrained block via Schur complement, evolve the reduced operator,
recover the constrained fields algebraically at the end. Both
formulas are on the slide.
