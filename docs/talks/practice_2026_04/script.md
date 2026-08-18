# Practice talk — spoken script

═══════════════════════════════════════════════════════════════
▶ Title slide (~15 s)
═══════════════════════════════════════════════════════════════

Hi everyone. I'm William, a master's student at the Cavendish.
Supervised by Will Barker, Mike Hobson, and Anthony Lasenby, I'm
working on whether extensions to our theory of gravity can amplify
a particular old graviton-photon conversion process. This is a
practice talk on my work the past few months on this project —
feedback on both content and delivery is very welcome.

▶ NEXT — slide 1: motivation

═══════════════════════════════════════════════════════════════
▶ Slide 1 — 01 motivation (~75 s)
═══════════════════════════════════════════════════════════════

The starting point is the Gertsenshtein effect: a graviton
propagating through a background magnetic field can, perhaps
surprisingly, convert into a photon.

Unfortunately, this effect is astrophysically useless. For a
magnetar with the strongest known magnetic fields, the conversion
probability is around ten-to-the-minus-ten. And for any
laboratory setup it's even more hopeless. At these scales, there
is no foreseeable detector that will see Gertsenshtein in any
astrophysical or laboratory setting... within standard-GR.

So our theoretical curiosity prompts us the question: is there
any extension of GR that can amplify this conversion to a
detectable level?

The approach in this project is unapologetically phenomenological
— instead of picking one specific extension of GR and defending
it, I'm systematically scanning the space of curvature, torsion,
and electromagnetic modifications of linearized GR and asking
which sectors, if any, amplify the conversion.

▶ NEXT — slide 2: the Gertsenshtein effect

═══════════════════════════════════════════════════════════════
▶ Slide 2 — 02 Gertsenshtein (~60 s)
═══════════════════════════════════════════════════════════════

I'm guessing you'll all first be wondering what the Gertsenshtein
effect actually is? In short: it's a mixing between the graviton
and the photon, mediated by the background magnetic field.

Although it might sound exotic, this isn't really new physics —
the coupling is already hiding in the Einstein–Maxwell action!
Once you turn on a background magnetic field and let the metric
fluctuate, those metric factors raising the field strength
indices produce a coupling between the graviton, the photon,
and the background field.

The equations are lengthy to reach, causing the surrounding
literature to be littered with mistakes, but their structure
is simple. It's a two-state mixing between the graviton and the
photon, with an off-diagonal coupling proportional to the
background magnetic field. The conversion becomes a standard
mixing process, a sin-squared oscillation with a beat frequency
set by that coupling.

▶ NEXT — slide 3: why it's so weak

═══════════════════════════════════════════════════════════════
▶ Slide 3 — 03 literature (~50 s)
═══════════════════════════════════════════════════════════════

But why is that coupling so weak? Fundamentally, it's the weakness
of gravity, with the conversion probability suppressed by the
square of the Planck mass.

Gertsenshtein himself in 1962 gave the formula and immediately
recognised it was hopelessly small. Boccaletti and collaborators
in 1970 did the full analytic solution for a localized B-field
region. Raffelt and Stodolsky in 1988 embedded the problem in
the axion-photon-mixing framework, which is the modern standard
for the effect.

We must remember that this formula is exact in linearized
Einstein–Maxwell — there are no approximations to wave away. Any
amplification has to come from genuine new physics.

▶ NEXT — slide 4: PGT and torsion

═══════════════════════════════════════════════════════════════
▶ Slide 4 — 04 pgt_a (~60 s)
═══════════════════════════════════════════════════════════════

Let me attempt to motivate the grounding of the theory we are
working with. Naturally, there exist two independent spacetime
symmetries: Lorentz invariance and translations. These are
fundamental, and we would expect that we could promote these
symmetries to local ones.

This is what a gauge theory does: it takes a global symmetry
and promotes it to a local one. So we have some new object which
describes how we do this at each point in spacetime, and we need
to add these compensating fields to ensure we don't mess up the
derivative terms.

It turns out that gauging the U(1) symmetry gives us precisely
the familiar field strength of QED. And the same construction
applied to gauging the Lorentz part of the Poincaré group gives
geometric curvature; applied to translations, it gives this
object called torsion.

As we know, GR keeps just the curvature, but at this level it
looks like there is good reason to keep both in our theory, and
with the Standard Model we see that gauge theories have found
great success in describing the other forces of nature.

As we do in GR, it's possible to assign a geometric intuition
to these objects. Parallel-transport a vector around a closed
loop and you get two mismatches: a rotational mismatch you know
as curvature, and a translational one where the loop literally
fails to close, which is torsion. There is also non-metricity,
describing the change in the length of that vector, which
neither GR nor PGT keeps.

▶ NEXT — slide 5: the parameter landscape

═══════════════════════════════════════════════════════════════
▶ Slide 5 — 05 pgt_b (~60 s)
═══════════════════════════════════════════════════════════════

Although we have a nice moral standpoint for PGT, we realize that
leaving the geometric perspective for a field theoretic viewpoint,
it becomes clear that there are a huge number of possible terms
that become available to construct our Lagrangian from.

There exists literature classifying the effects and phenomenology
of some of these terms, but the full landscape is largely
unexplored. We are free to add parity-odd terms, higher-derivative
terms, non-minimal couplings between fields, and so on. Not only
this but there ends up being many indices through which we can
contract over, so I had to write appropriate mathematica code to
systematically enumerate all the possible independent terms.

The table here already shows a vast number of terms, few of which
have been studied in the literature. And we don't even have to
stop there, with no reason a priori to exclude even higher order terms.

We will also consider more phenomenological extensions, such as
plasma effects, which are more like effective interactions than
fundamental modifications of gravity, but which are astrophysically
relevant.

▶ NEXT — slide 6: the question and the method

═══════════════════════════════════════════════════════════════
▶ Slide 6 — 06 question (~45 s)
═══════════════════════════════════════════════════════════════

The central question is therefore whether there exists within
this vast landscape any combination of these modifications that
can amplify the Gertsenshtein conversion. We are not looking for
a specific mechanism, but rather asking whether any mechanism
exists at all.

For every Lagrangian we test, we first derive the linearized
field equations symbolically — by hand they're intractable for
anything but the simplest case. We then route those equations
through our solver to compute the conversion probability and
compare it to the GR baseline, given they first pass stability
and causality checks.

This is what the bulk of my project has been about — building
the infrastructure to allow us to do this systematically and
efficiently.

▶ NEXT — slide 7: TIDAL overview

═══════════════════════════════════════════════════════════════
▶ Slide 7 — 07 tidal_overview (~60 s)
═══════════════════════════════════════════════════════════════

We start with a Lagrangian declaration, including the field
content, the nature of any perturbative expansion, and the
background. We then do the symbolic work to derive the equations
of motion, perturb them around the background, and put them in
a form suitable for numerical solution. Finally, we run the
numerical solver and do inference on the results.

I had to ensure that the framework was flexible enough to handle
the full range and structure of theories we want to test, and
that the symbolic and numerical layers were fully automated and
connected. Importantly for the inference layer, the numerical
solver has to have fast enough turnaround to allow us to run
thousands of simulations across the parameter space, and accurate
too.

▶ NEXT — slide 8: the modal solver

═══════════════════════════════════════════════════════════════
▶ Slide 8 — 08 tidal_numerical (~70 s)
═══════════════════════════════════════════════════════════════

The core numerical engine is the spectral solver tailored for the
linearized gauge theories the survey produces. The linearization
is key, allowing us to work in Fourier space and solve for each
mode independently, which is a huge speedup compared to a full
nonlinear PDE solver. I want to be able to handle any new terms
you want to add, and gauge theories generate many non-trivial
structures (constraint fields, near-degeneracies, etc.) that need
careful handling, so I had to make sure the solver was robust
enough to handle all of these.

Unfortunately though, the exact matrix solution is met with many
catastrophic numerical issues when applied to gauge theories,
which is what I had to fix to make it work for our survey. We are
not actually able to eigendecompose the operator, so we have to
compute the matrix exponential directly in a way that avoids the
ill-conditioned and hugely degenerate eigenvectors. I also had to
deal with the presence of algebraically constrained fields, which
enter the equations of motion but don't have their own explicit
dynamics.

The key benefit of the modal solver is that we don't actually end
up having to do any time-stepping at all, providing an exact
solution via a matrix exponential. That's a huge advantage in
terms of the machine precision accuracy of the time evolution,
which is crucial for the weak signals we're trying to measure.
The solver is also very fast, and the cost is independent of how
long I run the simulation, since it's a single matrix exponential
step that feeds directly into the inference layer.

▶ NEXT — slide 9: the inference framework

═══════════════════════════════════════════════════════════════
▶ Slide 9 — 09 inference (~75 s)
═══════════════════════════════════════════════════════════════

The question I was faced with when I started the project was:
how do I actually report results from this survey? I have a huge
parameter space, and for each point in that space there's a
conversion probability. How do we turn that into a statement
about whether amplification exists or not?

To this end, I was able to build a Bayesian inference layer on
top of the numerical solver, which allows us to turn the amplification
factor across the parameter space into a likelihood for nested
sampling. This gives us two key summary numbers: log Z, which
tells us the expected amplification for the theory, and the KL
divergence, which tells us whether any enhancement is broadly
distributed or concentrated within the parameter space.

We must keep to the linearized regime, taking the limit of
vanishingly small magnetic fields at which the amplification
becomes independent of the experimental setup and is purely a
characteristic property of the theory we wish to measure. In
precisely this regime we can also thankfully neglect any
backreaction on the geometry, and the problems that comes with
a uniform magnetic field filling all of spacetime.

▶ NEXT — slide 10: pre-campaign results

═══════════════════════════════════════════════════════════════
▶ Slide 10 — 11 results (~90 s)
═══════════════════════════════════════════════════════════════

To help characterize the direction of the full campaign, I ran
a few preliminary tests on some simple theories. The aim is to
not just run a single campaign with the most general theory, but
to also investigate the structure of theories of interest to the
community.

The first key model is the minimal PGT extension, which includes
the three torsion-mass invariants and a higher order R-squared
term. The Ricci-Cartan tensor promotes some torsion components
into propagating fields, but the equations of motion show a
complete structural decoupling between the new sectors and the
Gertsenshtein channel — the new torsion dynamics simply don't
enter the channel that an incoming graviton can excite. It was
clear that any amplification must come from non-minimal sectors
— meaning couplings that connect the fields in new ways.

Second was a more phenomenologically motivated theory, that of
pure Einstein–Maxwell plus an effective photon mass. The photon
mass detunes the graviton-photon resonance and suppresses the
conversion.

Another key area of interest in BSM phenomenology is the dark-photon
model, which is a hidden U(1) Proca vector kinetically mixed with
electromagnetism, and is a popular candidate for dark matter. Here,
I identify torsion's trace vector as this hidden field. In vacuum,
the relevant photon and dark photon sectors are exact eigenmodes
and no graviton-induced channel opens. So I added plasma to break
the degeneracy, but found the dark photon itself had negligible
effect on total conversion. Only the photon mass matters, and as
suppression.

▶ NEXT — slide 11: summary

═══════════════════════════════════════════════════════════════
▶ Slide 11 — 13 summary (~45 s)
═══════════════════════════════════════════════════════════════

Three things to take away.

First, standard Gertsenshtein is astrophysically useless — any
detectable signal needs new gravitational physics. Second, I'm
systematically scanning the quadratic curvature and torsion
landscape, with the inference cleanly separating amplification,
suppression, and fine-tuning. Third, the first half of the scan
is done — simple PGT extensions don't amplify, the non-minimal
Ricci-EM direction gives strong suppression, and the
higher-dimensional, parity-odd, and plasma corners are next.

The line on the slide is the one I want to land carefully: whether
or not amplification exists in this landscape, I'm building the
systematic map of where it could live.

Thanks — questions and delivery feedback are both welcome.
