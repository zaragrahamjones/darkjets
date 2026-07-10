===================================
INITIALIZATION:
===================================

darkinvisible is based on main507 in the PYTHIA examples, with the addition of a lund.dat output file for the jet analysis. The .cmnd file is used to set up the Pythia parameters for the simulation. The .cc file contains the main program that initializes Pythia, runs the event loop, and performs jet clustering using FastJet. The Makefiles are directly copied from the example folder too and Makefile.inc has been updated to include darkinvisible as a target for compilation.

plot_lund.py is a simple script to plot the lund.dat output file. It uses matplotlib to create a scatter plot of the log(1/Delta) vs log(kt) for the jets in the events.

===================================
COMPARISON WITH SM MODELS:
===================================

Within the .cmnd file, Zp mostly changes the rate and width, not a clean "amount of DM".

A clean SM comparison model would be something that gives the same _visible_ final state: a QCD jet recoiling against missing energy (neutrinos, anything else?). 
This could be achieved by:
- a Z' with the same mass and width, but decaying to SM particles instead of DM,
- pp-> Z -> QCD jets + invisible particles, and adjusting Z'=Z

PROBLEM: None of the produced DM particles decay back into SM particles, so final state does not actually contain dark jets, just missing energy. This is the scenario for a long lasting dark sector, but not what we are investigating.


===================================
# Pythia 8 Setup for Visible and Slightly Invisible Dark Jets

_Given to me by Mihoko, from ChatGPT_

This note summarizes a practical Pythia 8 setup for generating visible dark jets, and a simple extension in which dark hadrons can decay partly into invisible particles. The intended use is a phenomenological toy benchmark based on the Pythia 8 Hidden Valley module.

## 1. Basic idea

A visible dark jet can be modeled in Pythia 8 by enabling a Hidden Valley dark shower and forcing the resulting dark hadrons to decay into Standard Model particles.

The basic chain is:

```text
pp -> dark quark pair -> dark shower -> dark hadrons -> SM visible particles
```

For a mostly visible but slightly invisible model, one can allow a small branching fraction of the dark hadrons to decay into stable invisible particles:

```text
pi_d -> q qbar          with large branching ratio
pi_d -> chi chi         with small branching ratio
```

Here `chi` is a stable invisible particle introduced for phenomenological purposes.

## 2. Minimal visible dark jet card

The following is a simple Pythia 8 Hidden Valley card for prompt visible dark jets.

```text
! ------------------------------------------------------------
! pp -> dark quarks -> dark shower -> visible dark hadrons
! Pythia8 HiddenValley example card
! ------------------------------------------------------------

Beams:idA = 2212
Beams:idB = 2212
Beams:eCM = 13600.

! --- Standard QCD/QED settings
PartonLevel:ISR = on
PartonLevel:FSR = on
PartonLevel:MPI = on
HadronLevel:Hadronize = on

! --- Hidden Valley production
! Option A: internal production of dark quark pairs
HiddenValley:qqbar2DvDvbar = on
HiddenValley:gg2DvDvbar = on
HiddenValley:alphaOrder = 1

! --- Hidden Valley shower and hadronization
HiddenValley:FSR = on
HiddenValley:fragment = on
HiddenValley:Ngauge = 3
HiddenValley:nFlav = 1
HiddenValley:spinFv = 0

! Dark confinement scale
HiddenValley:Lambda = 10.0

! Shower cutoff
HiddenValley:pTminFSR = 1.0

! --- Dark quark mass
4900101:m0 = 20.0

! --- Dark meson masses
! pi_v
4900111:m0 = 10.0
! rho_v
4900113:m0 = 25.0

! --- Prompt visible decays
4900111:mayDecay = on
4900111:tau0 = 0.0
4900111:oneChannel = 1 1.0 91 1 -1

4900113:mayDecay = on
4900113:tau0 = 0.0
4900113:oneChannel = 1 1.0 91 1 -1
```

The key point is that the dark mesons should not be left stable if the goal is to generate visible dark jets. For example, `4900111` and `4900113` are allowed to decay promptly into Standard Model quarks.

## 3. Using an external hard process

For analysis work, it is often cleaner to generate the hard process externally, for example with MadGraph,

```text
pp -> Z' -> q_dark q_darkbar
```

and then let Pythia handle only the dark shower, fragmentation, and decays.

In that case, the Pythia card can start from an LHE file:

```text
Beams:frameType = 4
Beams:LHEF = events.lhe

HiddenValley:FSR = on
HiddenValley:fragment = on
HiddenValley:Ngauge = 3
HiddenValley:nFlav = 1
HiddenValley:spinFv = 0
HiddenValley:Lambda = 10.0
HiddenValley:pTminFSR = 1.0

4900101:m0 = 20.0

4900111:m0 = 10.0
4900111:mayDecay = on
4900111:tau0 = 0.0
4900111:oneChannel = 1 1.0 91 1 -1

4900113:m0 = 25.0
4900113:mayDecay = on
4900113:tau0 = 0.0
4900113:oneChannel = 1 1.0 91 1 -1
```

This separation is usually preferable when the mediator model, production cross section, or kinematics need to be controlled explicitly.

## 4. Meaning of the main parameters

| Setting | Meaning |
|---|---|
| `HiddenValley:FSR = on` | Enables final-state radiation in the dark sector. |
| `HiddenValley:fragment = on` | Enables dark-sector fragmentation / hadronization. |
| `HiddenValley:Ngauge = 3` | Chooses an SU(3)-like dark gauge group. |
| `HiddenValley:nFlav` | Number of dark quark flavors. |
| `HiddenValley:Lambda` | Dark confinement scale. |
| `HiddenValley:pTminFSR` | Dark shower cutoff scale. |
| `4900101:m0` | Dark quark mass. |
| `4900111:m0` | Dark pion mass. |
| `4900113:m0` | Dark rho mass. |
| `tau0 = 0.0` | Prompt decay. |
| `oneChannel` / `addChannel` | Defines decay channels and branching ratios. |

## 5. Adding a small invisible branching fraction

A simple way to move from a fully visible dark jet to a mostly visible / slightly invisible dark jet is to allow the dark pion to decay invisibly with a small branching ratio.

For example:

```text
pi_d -> q qbar      90%
pi_d -> chi chi     10%
```

The corresponding Pythia decay setup is:

```text
! ------------------------------------------------------------
! Mostly visible dark shower toy model
! pi_d -> SM quarks with 90%, invisible chi chi with 10%
! ------------------------------------------------------------

Beams:idA = 2212
Beams:idB = 2212
Beams:eCM = 13600.

! --- Hidden Valley production
HiddenValley:qqbar2DvDvbar = on
HiddenValley:gg2DvDvbar = on

! --- Hidden Valley shower and fragmentation
HiddenValley:FSR = on
HiddenValley:fragment = on
HiddenValley:Ngauge = 3
HiddenValley:nFlav = 2
HiddenValley:spinFv = 0
HiddenValley:Lambda = 10.0
HiddenValley:pTminFSR = 1.0

! --- Dark quarks
4900101:m0 = 20.0
4900102:m0 = 20.0

! --- Invisible stable particle chi
4900022:all = chi chi 2 0 0 1.0 0.0 0.0 0.0 0.0
4900022:mayDecay = off
4900022:isVisible = false

! --- Dark pion
4900111:m0 = 10.0
4900111:mayDecay = on
4900111:tau0 = 0.0

! 90% visible decay to light quarks
4900111:oneChannel = 1 0.90 91 1 -1

! 10% invisible decay to chi chi
4900111:addChannel = 1 0.10 100 4900022 4900022
```

Here `4900022` is introduced as a stable invisible particle. The line

```text
4900022:mayDecay = off
```

makes it stable in Pythia, and

```text
4900022:isVisible = false
```

marks it as invisible for event-record purposes.

## 6. Splitting visible decays among several quark flavors

Instead of forcing all visible decays into `d dbar`, one may split the visible branching ratio among light quarks:

```text
! pi_d visible 90%, invisible 10%

4900111:oneChannel = 1 0.30 91 1 -1
4900111:addChannel = 1 0.30 91 2 -2
4900111:addChannel = 1 0.30 91 3 -3
4900111:addChannel = 1 0.10 100 4900022 4900022
```

A more heavy-flavor-rich version could be:

```text
! b-rich visible dark pion

4900111:oneChannel = 1 0.20 91 1 -1
4900111:addChannel = 1 0.20 91 3 -3
4900111:addChannel = 1 0.50 91 5 -5
4900111:addChannel = 1 0.10 100 4900022 4900022
```

If `b bbar` decays are used, the dark pion mass should be safely above the `2 m_b` threshold. In practice, values such as

```text
4900111:m0 = 15.0
```

or larger are safer than a mass very close to threshold.

## 7. Interpretation of the invisible branching ratio

The invisible branching ratio specified in the decay table is a per-dark-hadron branching ratio. It is not exactly the same quantity as the event-level invisible energy fraction usually denoted by

```text
r_inv
```

A rough event-level quantity is

```text
r_inv ~ sum(E_invisible) / sum(E_dark_shower)
```

but the actual value depends on the dark hadron multiplicity, mass spectrum, boost, visible decay modes, neutrinos from Standard Model decays, and jet clustering.

Therefore, for a first phenomenological scan, it is useful to try

```text
BR(pi_d -> invisible) = 0.00, 0.05, 0.10, 0.20
```

and then measure the resulting event-level missing energy fraction in the generated sample.

## 8. Including dark rho decays

The dark rho can also be given visible and invisible decay channels:

```text
4900113:m0 = 25.0
4900113:mayDecay = on
4900113:tau0 = 0.0

4900113:oneChannel = 1 0.45 91 1 -1
4900113:addChannel = 1 0.45 91 3 -3
4900113:addChannel = 1 0.10 100 4900022 4900022
```

However, a more dark-shower-like setup may be to let the dark rho decay first into dark pions, and then let the dark pion decay visibly or invisibly:

```text
! rho_d -> pi_d pi_d
4900113:oneChannel = 1 1.0 100 4900111 4900111

! pi_d -> visible/invisible
4900111:oneChannel = 1 0.90 91 1 -1
4900111:addChannel = 1 0.10 100 4900022 4900022
```

This realizes a cascade in which a small invisible component appears at the end of the dark hadron decay chain.

## 9. Displaced visible or semi-visible decays

For prompt decays, one can use

```text
4900111:tau0 = 0.0
```

For displaced decays, set a nonzero lifetime and make sure Pythia is allowed to decay particles with that lifetime:

```text
ParticleDecays:limitTau0 = off
ParticleDecays:limitCylinder = on
ParticleDecays:xyMax = 30000
ParticleDecays:zMax = 30000

4900111:tau0 = 10.0
4900113:tau0 = 10.0
```

In common Pythia usage, `tau0` is effectively treated as `c tau` in mm. The precise treatment should be checked against the Pythia version and event-record output being used.

## 10. Suggested first benchmark

For a first mostly visible dark jet benchmark with a small invisible component, a compact setup is:

```text
HiddenValley:Ngauge = 3
HiddenValley:nFlav = 2
HiddenValley:Lambda = 10.0
HiddenValley:pTminFSR = 1.0

4900101:m0 = 20.0
4900102:m0 = 20.0

4900111:m0 = 15.0
4900111:mayDecay = on
4900111:tau0 = 0.0

4900022:all = chi chi 2 0 0 1.0 0.0 0.0 0.0 0.0
4900022:mayDecay = off
4900022:isVisible = false

4900111:oneChannel = 1 0.45 91 1 -1
4900111:addChannel = 1 0.45 91 3 -3
4900111:addChannel = 1 0.10 100 4900022 4900022
```

This corresponds to a dark pion with 90% visible decays and 10% invisible decays.

A useful first scan would be:

```text
BR_inv = 0.00, 0.05, 0.10, 0.20
m_pi_d = 10, 15, 20 GeV
Lambda_d = 5, 10, 20 GeV
```

with fixed dark quark and mediator masses, followed by measuring the generated distributions of jet mass, charged multiplicity, missing transverse momentum, and the alignment between the dark jet and missing momentum.

## 11. Caveats

The cards above are intended as toy-model benchmarks. For a publishable or reinterpretation-level study, the following points should be made explicit:

- the mediator model and production process,
- dark quark masses,
- dark confinement scale,
- dark pion and dark rho spectrum,
- visible and invisible branching ratios,
- particle lifetimes,
- whether the invisible fraction is defined at the decay-table level or measured at the event level,
- whether the analysis uses Pythia-internal production or an external LHE sample.

The most robust workflow is usually:

```text
MadGraph or another generator: hard process
Pythia 8: dark shower, dark fragmentation, and decays
Analysis code: event-level invisible fraction and jet observables
```
