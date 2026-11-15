---
title: NRTL
---

[TOC]

# NRTL

The Non-Random-Two-Liquid (NRTL) model presented by Renon and Prausnitz is an
activity-coefficient model for the excess Gibbs energy of liquid mixtures. It
uses local compositions to represent the non-randomness of the mixture.

Use NRTL when you have binary interaction parameters fitted for your components
and the property of interest (for example vapor-liquid or liquid-liquid data).
It is useful for strongly nonideal liquids, including systems that may demix,
but does not discover parameters from chemical names or molecular structure.
For structure-based estimates, see the [UNIFAC models](unifac_models/index.html).

Intuitively, the neighbors around a molecule need not have the same composition
as the bulk liquid: favorable contacts can be overrepresented. NRTL uses a
directional interaction measure \(\tau\) and a non-randomness parameter
\(\alpha\) to weight those local contacts. The result is an excess Gibbs energy,
relative to ideal liquid mixing, not a residual ideal-gas departure. Its
derivatives give \(\ln\gamma_i\); the dimensionless activity is
\(x_i\gamma_i\), not \(\gamma_i\) alone. See the
[shared introduction](index.html) for the distinction from fugacity coefficients.

# Theory

The model averages interaction energies around each component, using local
composition weights rather than just bulk mole fractions. With
\(n=\sum_i n_i\) and \(x_i=n_i/n\), its excess Gibbs energy is

\[
   G^E = nRT \sum_i x_i \frac{\sum_j x_j \tau_{ji} G_{ji}}{\sum_j x_j G_{ji}}
\]

with

\[
   \tau_{ij} = A_{ij} + \frac{B_{ij}}{T}, \qquad
   G_{ij} = \exp(-\alpha_{ij}\tau_{ij})
\]

where \(A_{ij}\), \(B_{ij}\), and \(\alpha_{ij}\) are the model parameters.
The negative sign and the product in \(\exp(-\alpha_{ij}\tau_{ij})\) are
essential; this is not \(\exp(-\alpha_{ij})\) or \(\exp(-\tau_{ij}/T)\).
In `yaeos`, the matrices `a`, `b`, and `c` store \(A_{ij}\), \(B_{ij}\), and
\(\alpha_{ij}\), respectively.

The activity coefficient of component \(i\) follows from

\[
   \ln \gamma_i =
   \frac{\sum_j x_j \tau_{ji} G_{ji}}{\sum_j x_j G_{ji}} +
   \sum_k \frac{x_k G_{ik}}{\sum_j x_j G_{jk}}
   \left(\tau_{ik} - \frac{\sum_j x_j \tau_{jk} G_{jk}}{\sum_j x_j G_{jk}}\right)
\]

# Parameters

The model is constructed from three \(N_c \times N_c\) matrices:

| Matrix | Meaning |
|:-------|:--------|
| `a` | Dimensionless \(A_{ij}\) interaction parameters |
| `b` | \(B_{ij}\) temperature-dependent interaction parameters [K] |
| `c` | Dimensionless \(\alpha_{ij}\) non-randomness parameters |

`T` is in kelvin. Store matrix element \((i,j)\) in `a(i,j)`, `b(i,j)`, and
`c(i,j)`; the column `(:,i)` enters the local average around component `i` in
the equation above. The two directions generally differ. Keep the same order
in `n` and every matrix, and normally set diagonal interaction coefficients to
zero. Non-randomness parameters are often symmetric, but must follow your fit.

Published NRTL tables do not all use this convention. If a source defines
\(\tau_{ij}=\Delta g_{ij}/(RT)\) with a constant energy difference, convert it
to `b(i,j) = Delta_g(i,j)/R` and set `a=0`, using a gas constant with matching
energy units. Do not enter J/mol directly as a kelvin-valued `b`, add another
minus sign without checking the source definition, or reuse UNIQUAC parameters.

# Usage

This complete program illustrates matrix orientation and the ideal-interaction
limit. The nonzero numbers are **illustrative**, not fitted parameters for any
named chemical pair.

```fortran
program nrtl_activity
use yaeos, only: pr, NRTL
implicit none

type(NRTL) :: model
real(pr) :: a(2, 2)
real(pr) :: b(2, 2)
real(pr) :: c(2, 2)
real(pr) :: n(2)
real(pr) :: T
real(pr) :: ln_gamma(2)
real(pr) :: ln_gamma_scaled(2)

a = 0.0_pr
b = 0.0_pr
b(1, 2) = 200.0_pr
b(2, 1) = -100.0_pr
c = 0.3_pr

model = NRTL(a, b, c)

n = [0.4_pr, 0.6_pr]
T = 300.0_pr
call model%ln_activity_coefficient(n, T, ln_gamma)
call model%ln_activity_coefficient(10.0_pr*n, T, ln_gamma_scaled)
print *, "ln(gamma), components 1 and 2 = ", ln_gamma
print *, "gamma = ", exp(ln_gamma)
print *, "Scaling difference, near zero = ", ln_gamma_scaled - ln_gamma

b = 0.0_pr
model = NRTL(a, b, c)
call model%ln_activity_coefficient(n, T, ln_gamma)
print *, "Zero-interaction ln(gamma), near zero = ", ln_gamma
end program nrtl_activity
```

`exp(ln_gamma)` is positive and may exceed 1; negative logarithms are normal.
Scaling all amounts should leave coefficients unchanged. When both `a` and `b`
vanish, this NRTL model becomes ideal (all coefficients are 1), regardless of
`c`. Fortran `reshape` fills columns first; explicit off-diagonal assignments
above avoid silently transposing directional parameters.

# Worked example: a liquid-liquid split

An activity model earns its keep when it predicts that a mixture separates. The
following program uses NRTL for a strongly nonideal binary and searches for a
two-liquid split with the library's stability and flash routines.

The parameters here are **illustrative**: they produce a split, but they are not
a fit to any named chemical pair.

```fortran
program nrtl_lle
use yaeos, only: pr, NRTL, EquilibriumState, flash, min_tpd
implicit none

integer, parameter :: nc = 2

type(NRTL) :: model
type(EquilibriumState) :: result
real(pr) :: a(nc, nc)
real(pr) :: b(nc, nc)
real(pr) :: c(nc, nc)
real(pr) :: n(nc)
real(pr) :: w(nc)
real(pr) :: T
real(pr) :: P
real(pr) :: mintpd
integer :: iters

a = 0.0_pr
b = 0.0_pr
c = 0.0_pr
a(1, 2) = 5.0_pr
a(2, 1) = 5.0_pr
c(1, 2) = 0.3_pr
c(2, 1) = 0.3_pr
model = NRTL(a, b, c)

n = [0.3_pr, 0.7_pr]
T = 300.0_pr
P = 1.0_pr

! Minimum of the reduced tangent-plane distance
call min_tpd(model, n, P=P, T=T, mintpd=mintpd, w=w)
print *, "mintpd = ", mintpd
print *, "trial composition = ", w

if (mintpd < 0.0_pr) then
   ! The trial composition gives the initial K-factors
   result = flash(model, n, T, k0=w/n, iters=iters)
   print *, "kind = ", result%kind
   print *, "x = ", result%x
   print *, "y = ", result%y
   print *, "beta = ", result%beta
else
   print *, "Feed is stable in this search."
end if
end program nrtl_lle
```

How to read this:

- `min_tpd` minimizes a tangent-plane distance measure at fixed `T` and `P`. A
  **negative** `mintpd` means the feed is unstable to a small perturbation, i.e.
  a second phase wants to appear. A positive value means the search found no
  instability; it is not proof that the feed is globally stable.
- `flash` for an activity model solves the isofugacity and mass-balance
  equations at a specified temperature. It **requires** the initial K-factors
  `k0`; without them the routine stops with an error rather than guessing.
  Use the trial composition from `min_tpd` for that purpose.
- `result%kind` reports the outcome. Only read `x`, `y`, and `beta` when it
  says `"split"`. The pressure and the two phase volumes are **not** set for an
  activity model and must not be read.
- `x` and `y` are the two liquid phases (for an activity model the
  "heavy"/"light" labels carry no density meaning, and `x` is not necessarily
  the larger phase), and `beta` is the mole fraction of the feed in the `y`
  phase. Check the material balance
  \(z_i=(1-\beta)x_i+\beta y_i\) yourself.

The same workflow is described in the
[phase-equilibrium guide](../phase_equilibrium/index.html); the
[flash page](../phase_equilibrium/flash.html) covers the activity-model branch in
detail.

# Notes

NRTL is a `GeModel`: it provides activity coefficients and excess thermal
properties, but no volumetric properties, and the flash routines solve only
the isofugacity and mass-balance equations.

Values such as 0.2 or 0.3 are common choices of \(\alpha_{ij}\). They are not
universal constants, and they are not an instruction to replace the value used
in a parameter fit: the non-randomness parameter is part of the fit, and
changing it changes the model's predictions.

The three matrices need not be symmetric, but the diagonal entries are usually
zero. The diagonal represents a component's interaction with itself, which
carries no information about unlike pairs.

Ordinary NRTL and [NRTL-HV](nrtl-hv.html) share the same mathematical form
when all covolumes are equal: with a unit covolume for every component,
NRTL-HV reproduces NRTL exactly when `gji0 = b_matrix * R` and
`gjiT = a_matrix * R`. Use this as a numerical cross-check when converting a
parameter set, not as a claim that the two models are interchangeable for
unequal covolumes.

## Common mistakes and next steps

Do not extrapolate a good room-temperature binary fit across a large
temperature range without checks. The temperature dependence of NRTL is
strictly linear in \(1/T\), so it cannot represent a change in interaction
enthalpy over the range, and the prediction drifts silently. Large values of
`a+b/T` can also overflow the exponential, and an apparently smooth curve can
still predict the wrong phase behavior. NRTL does not include ionic
speciation, chemical reaction, or density physics by itself.

A good sequence of work is: compare activity coefficients and excess enthalpies
against your parameter source using the
[shared property methods](index.html); then use
[stability and phase-equilibrium calculations](../phase_equilibrium/index.html)
to determine whether a composition is stable, remembering that a single call
to an activity method is not a flash calculation; and finally compare the
result with data. For vapor-liquid equilibrium you also need the liquid
standard-state and vapor fugacity information, not just activity coefficients.
For a cubic EOS with Huron-Vidal mixing, read [NRTL-HV](nrtl-hv.html) rather
than assuming the two parameterizations are interchangeable.

# References

1. H. Renon and J. M. Prausnitz, "Local compositions in thermodynamic excess
   functions for liquid mixtures," AIChE Journal, vol. 14, no. 1, pp. 135-144,
   1968, doi: 10.1002/aic.690140124.
