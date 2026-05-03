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

I'm guessing the first question you'll all have is: what actually is the Gertsenshtein effect? The short answer is: it's a mixing between the graviton and the photon, mediated by the background magnetic field.

Though perhaps sounding exotic, this is not a new interaction vertex — it's already present in the Einstein–Maxwell action! If we recall that raising indices requires a factor of the inverse metric, we can see the action contains a sort-of three-field vertex coupling between the graviton, the photon, and the background magnetic field.

The equations are lengthy to reach, causing the surrounding literature to be littered with mistakes, but their structure is simple. It's a two-state mixing between the graviton and the photon, with an off-diagonal coupling proportional to the background magnetic field. The conversion becomes a standard mixing process, with a beat frequency set by that coupling.

▶ NEXT — slide 3: why it's so weak

═══════════════════════════════════════════════════════════════
▶ Slide 3 — 03 literature (~50 s)
═══════════════════════════════════════════════════════════════

But why is the conversion so weak? Its fundamental origin is the weakness of gravity, suppressed by the huge Planck scale. It's not exactly looking good for the gertsenshtein effect when the conversion probability is already suppressed by the Planck mass squared.

Gertsenshtein himself in 1962
gave the formula and immediately recognised it was hopelessly small.
Boccaletti and collaborators in 1970 did the full analytic solution for
a localised B-field region. Raffelt and Stodolsky in 1988 embedded the
problem in the axion-photon-mixing framework, which is the modern standard for the effect.

We must remember that this formula is exact in linearised Einstein–Maxwell, there are no approximations to wave away. Any amplification has to come from genuine new physics.

▶ NEXT — slide 4: PGT and torsion

═══════════════════════════════════════════════════════════════
▶ Slide 4 — 04 pgt_a (~60 s)
═══════════════════════════════════════════════════════════════

Let me attempt to motivate the grounding of the theory we are working with. Naturally, there exists in spacetime two potentially independent symmetries: Lorentz invariance and translations. These are fundamental, and we would expect that we could promote these symmetries to local ones.

This is what a gauge theory does: it takes a global symmetry and promotes it to a local one. So we have some new object which describes how we do this at each point in spacetime, and we need to add these compensating fields to ensure we dont mess up the derivative terms.

It turns out that gauging the U(1) symmetry gives us precisely the familiar field strength of QED. It turns out that the field strength arising from gauging the Lorentz part of the Poincaré group gives us exactly the geometric curvature, while gauging translations gives us this new object called torsion.

As we know, GR keeps just the curvature, but at this level it looks like there is good reason to keep both in our theory, and we already know from the Standard Model that gauge theories have found great success in describing the other forces of nature.

As we do in GR, its possible to assign a geometric intuition to these objects.
Parallel-transport a vector around a closed loop and you get two mismatches: a rotational
mismatch you know as curvature, and a translational one where the loop
literally fails to close, which is torsion. There is also non-metricity, describing the change in the length of that vector, which neither GR nor PGT keeps.

▶ NEXT — slide 5: the parameter landscape

═══════════════════════════════════════════════════════════════
▶ Slide 5 — 05 pgt_b (~60 s)
═══════════════════════════════════════════════════════════════

Although we have a nice moral standpoint for PGT, we realise that leaving the geometric perspective for a field theoretic viewpoint, it becomes clear that there are a huge number of possible terms that become available to construct our Lagrangian from.

There exists literature classifying the effects and phenomenology of some of these terms, but the full landscape is largely unexplored. We are free to add parity odd-terms, higher-derivative terms, non-minimal couplings between fields, and so on. Not only this but there ends up being many indices through which we can contract over, so i had to write appropriate mathematica code to systematically enumerate all the possible independent terms.

The table here already shows a vast number of terms, few of which have been studied in the literature. And we don't even have to stop it there, with no reason a priori to exclude even higher order terms.

We will also consider more phenomenological extensions, such as plasma effects, which more like effective interactions than fundamental modifications of gravity, but which are astrophysical relevant.

▶ NEXT — slide 6: the question and the method

═══════════════════════════════════════════════════════════════
▶ Slide 6 — 06 question (~45 s)
═══════════════════════════════════════════════════════════════

The central question is therefore whether there exists within this vast landscape any combination of these modifications that can amplify the Gertsenshtein conversion. We are not looking for a specific mechanism, but rather asking whether any mechanism exists at all.

For every Lagrangian we test, we first derive
the linearised field equations symbolically — by hand they're
intractable for anything but the simplest case. We then route those equations
through our solver to compute the conversion probability and compare it to the GR baseline.

This is what the bulk of my project has been about — building the infrastructure to do allow us to do this systematically and efficiently.

▶ NEXT — slide 7: TIDAL overview

═══════════════════════════════════════════════════════════════
▶ Slide 7 — 07 tidal_overview (~60 s)
═══════════════════════════════════════════════════════════════

We start with a Lagrangian declaration, including the field content, the nature of any perturbative expansion, and the background. We then do the symbolic work to derive the equations of motion, perturb them around the background, and put them in a form suitable for numerical solution. Finally, we run the numerical solver and do inference on the results.

I had to ensure that the framework was flexible enough to handle the full range and structure of theories we want to test, and that the symbolic and numerical layers were fully automated and connected. Importantly for the inference layer, the numerical solver has to fast enough turnaround to allow us to run thousands of simulations across the parameter space, and accurate too.

▶ NEXT — slide 8: the modal solver

═══════════════════════════════════════════════════════════════
▶ Slide 8 — 08 tidal_numerical (~70 s)
═══════════════════════════════════════════════════════════════

The core numerical engine is the spectral solver tailored for the
linearised gauge theories the survey produces. The linearisation is key, allowing us to work in Fourier space and solve for each mode independently, which is a huge speedup compared to a full nonlinear PDE solver.
We want a framework that will be able to handle any new terms you want to
add to it, and in general there are a lot of different structures and types of solutions that can arise, so we had to make sure the solver was flexible enough to handle all of these.

The key benefit of the modal solver is that we don't actually end up having
to do any time-stepping at all, provided an exact solution via a matrix exponential. That's a huge advantage in terms of the machine precision accuracy of the time evolution, which is crucial for the weak signals we're trying to measure. The solver is also very fast, and constant order in the end time, since it does only a single time step.

Unfortunately though, the exact matrix solution is met with many catastrophic numerical issues when applied to gauge theories, which is what we had to fix to make it work for our survey. We are not actually able to eigendecompose the operator, so we have to compute the matrix exponential directly in a way that avoids the ill-conditioned and hugely degenerate eigenvectors. We also have to deal with the presence of algebraically constrained fields, which enter the equations of motion but don't have their own explicit dynamics.

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
