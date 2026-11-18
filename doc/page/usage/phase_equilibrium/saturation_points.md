---
title: Saturation Points
---

[TOC]

# Saturation points

A saturation point is a state where a feed of fixed composition is in
equilibrium with an incipient phase. The two most common cases are:

- **Bubble point**: the incipient phase is vapor and the feed is liquid.
- **Dew point**: the incipient phase is liquid and the feed is vapor.
- **Liquid-liquid point**: both phases are liquids.

In `yaeos`, saturation points are calculated with
[[saturation_pressure]] and [[saturation_temperature]]. Both return an
[[EquilibriumState]].

“Incipient” means the limit as the amount of the new phase approaches zero.
The new phase can have a very different composition despite its negligible
amount. Saturation fixes the bulk-phase composition and **either** temperature
or pressure, then solves the other condition and the incipient composition.
To find finite amounts of both phases at specified T and P, use a
[flash](flash.html). These saturation interfaces require an `ArModel`.

# Theory

A saturation point is the limit of a two-phase region at fixed composition. For
a bubble point, the liquid phase has the feed composition and the vapor phase is
incipient:

$$
   z_i = x_i, \qquad \sum_i x_i K_i = 1
$$

For a dew point, the vapor phase has the feed composition and the liquid phase
is incipient:

$$
   z_i = y_i, \qquad \sum_i \frac{y_i}{K_i} = 1
$$

Here \(x_i\) and \(y_i\) denote liquid and vapor mole fractions and
\(K_i=y_i/x_i\). The bubble equation normalizes the emerging vapor,
`y=x*K`; the dew equation normalizes the emerging liquid, `x=y/K`.

In both cases the K-factors are updated from the fugacity coefficients, and the
pressure or temperature is adjusted until the corresponding sum condition is
satisfied. A liquid-liquid point follows the same structure with both phases
evaluated on the liquid root, so the sum condition normalizes the incipient
liquid composition.

# `saturation_pressure`

```fortran
type(EquilibriumState) function saturation_pressure( &
   model, n, T, kind, p0, y0, max_iters &
   )
   class(ArModel), target, intent(in) :: model
   real(pr), intent(in) :: n(:)          !! Composition [moles or mole fractions]
   real(pr), intent(in) :: T             !! Temperature [K]
   character(len=*), intent(in) :: kind  !! "bubble", "dew", or "liquid-liquid"
   real(pr), optional, intent(in) :: p0  !! Initial pressure [bar]
   real(pr), optional, intent(in) :: y0(:) !! Initial incipient composition
   integer, optional, intent(in) :: max_iters
end function saturation_pressure
```

The calculation first solves a scalar normalization equation. Internally its
updated ratio \(q_i\) is **incipient/bulk**, so the common equation is:

$$
    f(P) = \sum_i z_i q_i(P) - 1 = 0
$$

For bubbles \(q=K\), but for dews \(q=1/K\). Do not transfer this internal
ratio convention into flash `k0`, which is always `y/x`. The scalar loop
updates this ratio from bulk/incipient fugacity coefficients. A `y0` argument
is always an incipient **composition**, not a K-factor vector.

A one-dimensional Newton step updates pressure. If the scalar iteration
reaches its limit, the
routine switches to a full Newton method that solves for the K-factors and
\(\ln P\) simultaneously.

```fortran
program bubble_pressure
   use yaeos, only: pr, ArModel, EquilibriumState, PengRobinson76
   use yaeos, only: saturation_pressure
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriumState) :: sat
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: n(2)
   real(pr) :: T
   real(pr) :: K(2)

   tc = [190.564_pr, 425.12_pr]
   pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]

   model = PengRobinson76(tc, pc, w)

   n = [2.5_pr, 6.7_pr]
   T = 150.0_pr

   sat = saturation_pressure(model, n, T, kind="bubble")
   if (sat%kind /= "bubble") error stop "No accepted bubble point"
   if (.not. ieee_is_finite(sat%P)) error stop "Nonfinite saturation pressure"
   if (sat%P <= 0.0_pr) error stop "Invalid saturation pressure"

   print *, sat%kind, sat%T, sat%P
   print *, "x = ", sat%x
   print *, "y = ", sat%y
   print *, "Composition sums: ", sum(sat%x), sum(sat%y)
   K = sat%y/sat%x
   print *, "Bubble sum x*K (should approach one): ", sum(sat%x*K)
end program bubble_pressure
```

This sum verifies the convention but is not an independent equilibrium test.
Also evaluate the log-fugacity mismatch as in the [flash example](flash.html).
Check finite, positive conditions, distinct volumes, and normalized
compositions before accepting the result.

# `saturation_temperature`

```fortran
type(EquilibriumState) function saturation_temperature( &
   model, n, P, kind, t0, y0, max_iters &
   )
   class(ArModel), target, intent(in) :: model
   real(pr), intent(in) :: n(:)          !! Composition [moles or mole fractions]
   real(pr), intent(in) :: P             !! Pressure [bar]
   character(len=*), intent(in) :: kind  !! "bubble", "dew", or "liquid-liquid"
   real(pr), optional, intent(in) :: t0  !! Initial temperature [K]
   real(pr), optional, intent(in) :: y0(:) !! Initial incipient composition
   integer, optional, intent(in) :: max_iters
end function saturation_temperature
```

The same two-stage strategy is used, but now with respect to temperature.
Call fragments using a declared model and result:

```fortran
sat = saturation_temperature(model, n, P=15.0_pr, kind="bubble")
sat = saturation_temperature(model, n, P=15.0_pr, kind="dew", t0=330.0_pr)
```

The optional initial guesses `p0`, `t0`, and `y0` are useful for a sequence of
points. For bubbles use `y0=previous%y`; for dews use `y0=previous%x`, because
the argument means **incipient** composition in either case. Reuse only checked
points. Initialization can help branch tracking but does not guarantee it.

# Choosing the kind

| `kind` | Feed phase | Incipient phase | Typical use |
|:-------|:-----------|:----------------|:------------|
| `"bubble"` | Liquid | Vapor | Boiling point at fixed composition |
| `"dew"` | Vapor | Liquid | Condensation point at fixed composition |
| `"liquid-liquid"` | Liquid | Liquid | Immiscible liquid splits |

The incipient phase root is selected internally. For `"bubble"` the main phase
uses the liquid root and the incipient phase the vapor root; for `"dew"` the
roles are reversed.

# Results

The returned [[EquilibriumState]] contains the same fields used by flash
calculations:

| Field | Meaning |
|:------|:--------|
| `kind` | The requested kind, or `"failed"` |
| `x` | Heavy-phase composition: the feed for a bubble point, the incipient liquid for a dew point |
| `y` | Light-phase composition: the incipient vapor for a bubble point, the feed for a dew point |
| `Vx`, `Vy` | Phase molar volumes [L/mol] |
| `T`, `P` | Saturation temperature [K] and pressure [bar] |
| `beta` | Phase fraction (0 for bubble, 1 for dew) |
| `iters` | Iterations used |

A detected failure has `kind="failed"`. Current checks include NaN values,
an iteration-limit test, and a phase-volume difference below `1e-8`. They are
not exhaustive residual or global-stability tests. Although `max_iters` is
accepted, the current fallback and final limit test use the module's
`max_iterations` rather than that argument. Do not rely on it to bound work.

# Pure-component saturation lines

For a pure component, [[Psat_pure]] calculates the saturation pressure at a
given temperature:

The following are fragments for an initialized model. Add imports and
declarations before executable statements, and choose a subcritical temperature
for the selected component (300 K is not subcritical for methane):

```fortran
real(pr) :: Psat
real(pr) :: Vl
real(pr) :: Vv
logical :: converged

Psat = model%Psat_pure(1, T=150.0_pr, Vl=Vl, Vv=Vv, converged=converged)
if (.not. converged) error stop "Pure saturation did not converge"
```

The method uses a secant/Newton formulation on the isofugacity condition and
returns the liquid and vapor volumes. It has no critical-temperature guard: for
methane at 300 K, above its critical temperature, it returned a pressure of
about 34.8 bar with identical liquid and vapor volumes and still reported
`converged=.true.`. Always compare the temperature with the component's
critical temperature and check that the two volumes differ.

To calculate a complete pure-component saturation line, use
[[pure_saturation_line]]:

```fortran
use yaeos, only: pure_saturation_line, PurePsat

type(PurePsat) :: psat

psat = pure_saturation_line(model, component=1, minP=0.1_pr, minT=100.0_pr)

! Interpolate the saturation temperature at a given pressure
print *, psat%get_T(P=1.0_pr)

! Interpolate the saturation pressure at a given temperature
print *, psat%get_P(T=150.0_pr)
```

The [[PurePsat]] type stores the saturation temperatures, pressures, and
liquid/vapor volumes, together with interpolators for `get_T` and `get_P`.
The line is traced from the pure-component critical point down to `minP` or
`minT` with the continuation method. Both lower limits are required, and the
trace stops at whichever limit is reached first. For methane with
`minP=0.1 bar` and `minT=100 K` the line covered 99.8 K to 190.5 K, that is,
up to the critical point. Inspect the stored `P` and `T` ranges before
interpolating; extrapolation is not evidence that liquid-vapor coexistence
exists.

# Relation with phase envelopes

A saturation point is the starting point for a phase-envelope calculation. The
envelope routines in [Phase envelopes](envelopes.html) take a converged
`EquilibriumState` and continue it along the phase boundary:

Fragment (add imports and declarations before executable statements):

```fortran
use yaeos, only: pt_envelope_2ph, PTEnvel2

type(PTEnvel2) :: envelope

sat = saturation_temperature(model, n, P=1.0_pr, kind="dew", t0=150.0_pr)
if (sat%kind == "dew") envelope = pt_envelope_2ph(model, n/sum(n), sat)
```

This is usually more robust than calculating every saturation point
independently, especially near critical points.

# Critical points

Critical points and critical lines are calculated with
[[critical_point]] and [[critical_line]]. They are described in the
[critical points](critical.html) page. The `spec_CP` object provides the
allowed specifications for a critical-point calculation.

# Practical notes

- Use `z = n/sum(n)` when working with mole fractions; the routines normalize
  the feed internally, and the compositions returned are normalized.
- Prefer a reduced model containing only present components for feeds with
  exact zeros. Some paths divide by feed fractions or take logarithms before
  masking zeros; blanket zero-component support is not assured.
- Nearly identical volumes can trigger failure near criticality. Do not treat
  an identical-phase solution as a physical bubble/dew point; investigate with
  the [critical solver](critical.html).
- When scanning a phase boundary, reuse `y0` and the previous pressure or
  temperature as initial guesses.
- See `example/basics/2_saturation_points.f90` and
  `example/basics/5_pure_psat.f90` for complete runnable examples.
