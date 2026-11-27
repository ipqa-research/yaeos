---
title: Phase Equilibrium
ordered_subpage: flash.md
                 saturation_points.md
                 envelopes.md
                 stability.md
                 critical.md
                 multiphase.md
---

Imagine cooling a closed bottle containing methane and butane. Initially it
may contain only gas; after cooling, a liquid can appear. Both substances can
be present in **both** phases, but in different proportions.

- A **component** is a chemical substance in the model (methane, butane, ...).
- A **phase** is a homogeneous portion of the system (a vapor, a liquid, or a
  second immiscible liquid). Two components do not imply two phases.
- The **feed composition** describes the entire bottle, not either phase
  separately. A phase fraction tells you how much material is in a phase,
  whereas a phase composition tells you what that phase contains.

These routines predict equilibrium according to the selected model. They do
not describe how fast droplets form, heat transfer, or the time required to
reach equilibrium.

In `yaeos`, most phase-equilibrium procedures return the
[[EquilibriumState]] type, which holds the phase compositions, molar volumes,
pressure, temperature, phase fraction, and iteration information. Multiphase
calculations return [[MPEquilibriumState]].

# Thermodynamic background

Use mole fractions throughout: `z = n/sum(n)` converts component mole amounts
`n` to a normalized feed `z`. Keep the component order identical to the model's
critical constants and interaction parameters. Temperatures are in K,
pressures in bar, and phase molar volumes in L/mol.

For a mixture at fixed temperature and pressure, phase equilibrium requires
the equality of the chemical potentials, or equivalently the fugacities, of
each component in every phase:

$$
   x_i^L \hat{\phi}_i^L = x_i^V \hat{\phi}_i^V \qquad (i = 1, \dots, N_c)
$$

together with the overall mass balance

$$
   z_i = (1-\beta)\,x_i + \beta\,y_i \qquad (i = 1, \dots, N_c)
$$

where \(z_i\) is the feed mole fraction, \(x_i\) and \(y_i\) are the heavy- and
light-phase mole fractions, and \(\beta\) is the light-phase mole fraction.
Eliminating \(x_i\) and \(y_i\) with the K-factors \(K_i = y_i/x_i\) gives the
Rachford-Rice equation:

$$
   \sum_{i=1}^{N_c} \frac{z_i\,(K_i-1)}{1+\beta\,(K_i-1)} = 0
$$

Here \(N_c\) is the number of components, \(i\) selects one component, and
\(\hat\phi\) is its dimensionless fugacity coefficient, a measure of
non-ideal behavior. \(L\) and \(V\) label liquid and vapor. Equal fugacity
means there is no thermodynamic driving force to transfer a component between
the phases. \(K_i>1\) means component \(i\) is enriched in the `y` phase,
not that it is absent from the `x` phase.

The fraction \(\beta\) is on a **mole basis**, not a mass or volume basis.
For a feed of \(N\) moles, the phases contain \((1-\beta)N\) and
\(\beta N\) moles. Their volumes are therefore \((1-\beta)NV_x\) and
\(\beta NV_y\). A small vapor mole fraction can occupy a large volume.

For the `GeModel` flash and stability routines discussed here, activity
coefficients replace fugacity coefficients to describe **liquid-liquid**
equilibrium. These routines do not automatically add vapor pressures or a
vapor-phase model to create a vapor-liquid calculation.

[TOC]

# Calculation types

The main phase-equilibrium calculations are:

| Calculation | Specification | Routine |
|:------------|:--------------|:--------|
| Two-phase flash | \(zPT \rightarrow x, y, \beta\) | [[flash]], [[flash_second_order]] |
| Two-phase flash | \(zVT \rightarrow x, y, \beta, P\) | [[flash]] |
| Bubble/dew pressure | Feed composition and \(T\); solve \(P\) and incipient composition | [[saturation_pressure]] |
| Bubble/dew temperature | Feed composition and \(P\); solve \(T\) and incipient composition | [[saturation_temperature]] |
| Phase envelope | Continuation from a saturation point | [[pt_envelope_2ph]] |
| Stability analysis | Tangent-plane distance | [[min_tpd]] |
| Critical point | Critical conditions | [[critical_point]] |
| Multiphase flash | \(zPT \rightarrow\) several phases | [[pt_mp_flash]] |

The detailed pages are organized as follows:

- [Flash calculations](flash.html) describes the \(zPT\) and \(zVT\) flashes,
  initial K-factors, second-order flash, and flashes with activity-coefficient
  models.
- [Saturation points](saturation_points.html) describes bubble, dew, and
  liquid-liquid points, plus pure-component saturation lines.
- [Phase envelopes](envelopes.html) describes the continuation method for
  two-phase envelopes and isopleths.
- [Stability analysis](stability.html) describes tangent-plane distance
  calculations and how they are used to initialize flashes.
- [Critical points](critical.html) describes critical-point and critical-line
  calculations.
- [Multiphase flash](multiphase.html) describes the calculation of systems with
  more than two equilibrium phases.

## Which calculation should I choose?

- **“How much liquid and vapor do I have at this T and P?”** Use a flash.
- **“At what pressure does the first bubble/drop appear at this T?”** Use
  saturation pressure. Fix pressure instead to find a saturation temperature.
  At a boundary the new phase has zero limiting amount; saturation does not
  determine a finite split inside the two-phase region.
- **“Where is the boundary over a range of T and P?”** Trace an envelope from
  one checked saturation point. It follows a branch, not every possible
  boundary of a complicated mixture.
- **“Should this homogeneous feed split at all?”** Start with stability
  analysis. A failed two-phase flash is not the answer to that question.
- **“Where do the two phases become identical?”** Use critical calculations.
- **“Could two liquids coexist with vapor?”** Read the multiphase page,
  including the current wrapper's implementation limitations.

## A volume root is not a complete equilibrium result

An equation of state can have several volumes at the **same composition**,
temperature, and pressure. `root_type="liquid"` or `"vapor"` selects a root;
`"stable"` selects the preferred homogeneous root at that composition. This
does not test whether a lower-Gibbs-energy state exists at **different phase
compositions**.

Likewise, equality of fugacities and material balance give a candidate
equilibrium. A local solver can converge to identical phases, a metastable
split, or an equilibrium that is unstable to a third phase. Global equilibrium
requires the lowest appropriate thermodynamic potential among all admissible
states. Stability searches help test this but finite numerical searches are
not a proof that no other state exists.

# The `EquilibriumState` type

A two-phase result contains:

| Field | Meaning |
|:------|:--------|
| `kind` | Routine-dependent: `"split"`, `"bubble"`, `"dew"`, `"liquid-liquid"`, `"critical"`, or `"failed"` |
| `x` | Heavy-phase mole fractions |
| `y` | Light-phase mole fractions |
| `Vx`, `Vy` | Phase molar volumes [L/mol] |
| `T`, `P` | Temperature [K] and pressure [bar] |
| `beta` | Mole fraction of the light phase |
| `iters` | Iterations used |

# General workflow

A typical phase-equilibrium workflow is:

1. Build a model with the desired equation of state and mixing rule.
2. Check the stability of the feed with [[min_tpd]] when the phase behaviour
   is uncertain.
3. Calculate a flash at the conditions of interest.
4. Use the resulting compositions as initial guesses for nearby points.
5. For a complete phase diagram, start from a saturation point and continue it
   with the envelope routines.

Start with the complete [PT flash example](flash.html). Its material-balance
and fugacity checks are useful even when the solver reports success.

# Practical notes

- Check the routine's result contract, finite positive conditions, normalized
  nonnegative compositions, phase fractions, and component balance. Failure
  handling is not uniform; critical and multiphase routines need extra checks.
- Use initial guesses (`k0`, `p0`, `t0`, `y0`) from previous calculations when
  scanning a region. Continuation is more robust than independent calculations.
- Near critical points, compositions and volumes approach each other and phase
  fractions become poorly determined. Try continuation, smaller steps, and
  alternative initial guesses; changing solver does not remove this sensitivity.
- For liquid-liquid equilibria, perform a stability analysis first and consider
  forcing the light phase to the liquid root in [[flash_second_order]].
- Activity-coefficient models ([[GeModel]]) require an initial K-factor
  estimate for `flash`; [[min_tpd]] can provide one. `flash_second_order`
  initializes these models with stability analysis when `k0` is absent.
- Validate model parameters against data in the temperature, pressure, and
  composition region of interest. A converged cubic-model water/hydrocarbon
  split, for example, is not by itself a validated prediction of that system.
