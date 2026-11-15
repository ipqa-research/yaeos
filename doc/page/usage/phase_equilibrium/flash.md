---
title: Flash calculations
---

[TOC]

# Flash calculations

Flash calculations determine how a feed of known composition splits into two or
more equilibrium phases at specified conditions. They are one of the most common
phase-equilibrium calculations in process and reservoir modelling.

In `yaeos`, flash routines return the [[EquilibriumState]] derived type, which
contains the phase compositions, molar volumes, pressure, temperature, phase
fraction, and the number of iterations used.

The two main specifications are:

- \(zPT \rightarrow x, y, \beta\): `flash` with `p_spec` and `T`.
- \(zVT \rightarrow x, y, \beta, P\): `flash` with `v_spec` and `T`.

# Theory

A two-phase flash solves the isofugacity condition and the mass balance
simultaneously. With K-factors \(K_i = y_i/x_i\), the mass balance gives the
Rachford-Rice equation:

$$
   \sum_{i=1}^{N_c} \frac{z_i\,(K_i-1)}{1+\beta\,(K_i-1)} = 0
$$

which is solved for the light-phase fraction \(\beta\). The K-factors are updated
from the fugacity coefficients:

$$
   K_i = \frac{\hat{\phi}_i^L}{\hat{\phi}_i^V}
$$

Here `z` is the normalized feed composition, `x` and `y` are the heavy- and
light-phase compositions, and `beta` is the mole fraction of feed in the `y`
phase. Once `K` and `beta` are known,

$$
 x_i=\frac{z_i}{1+\beta(K_i-1)}, \qquad y_i=K_i x_i.
$$

Both phase compositions must sum to one. The balance
`z = (1-beta)*x + beta*y` must hold component by component, not just in total.

For an `ArModel`, fugacity coefficients are evaluated from the equation of
state. For a `GeModel`, these flash routines use activity coefficients in
both phases: this is a liquid-liquid calculation, not a vapor-liquid flash
with automatically supplied pure-component vapor pressures.

# `flash`

The [[flash]] function uses successive substitution. Its interface is:

```fortran
type(EquilibriumState) function flash(model, z, T, v_spec, p_spec, k0, iters)
   class(BaseModel), intent(in) :: model
   real(pr), intent(in) :: z(:)
   real(pr), intent(in) :: T
   real(pr), optional, intent(in) :: v_spec
   real(pr), optional, intent(in) :: p_spec
   real(pr), optional, intent(in) :: k0(:)
   integer, optional, intent(out) :: iters
end function flash
```

For an `ArModel`, provide **exactly one** of `p_spec` and `v_spec`. Pass
normalized mole fractions `z`; unlike saturation routines, this routine does
not normalize the feed for you. For a `GeModel`, omit both specifications and
provide `k0`. The result is a
`type(EquilibriumState)` with:

| Field | Meaning |
|:------|:--------|
| `kind` | `"split"` for an accepted split, `"failed"` for a trivial solution or iteration failure |
| `x` | Heavy-phase mole fractions |
| `y` | Light-phase mole fractions |
| `Vx` | Heavy-phase molar volume [L/mol] |
| `Vy` | Light-phase molar volume [L/mol] |
| `T` | Temperature [K] |
| `P` | Pressure [bar] |
| `beta` | Mole fraction of the light phase |
| `iters` | Iterations used |

## PT flash

```fortran
program pt_flash
   use yaeos, only: pr, ArModel, EquilibriumState, PengRobinson76, flash
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriumState) :: result
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: z(2)
   real(pr) :: T
   real(pr) :: P
   real(pr) :: balance(2)
   real(pr) :: lnphi_x(2)
   real(pr) :: lnphi_y(2)
   integer :: iters

   tc = [190.564_pr, 425.12_pr]  ! Methane / n-butane
   pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]

   model = PengRobinson76(tc, pc, w)

   z = [0.4_pr, 0.6_pr]
   T = 294.0_pr
   P = 60.0_pr

   ! Supply iters explicitly; see the implementation note below.
   result = flash(model, z, T, p_spec=P, iters=iters)
   if (result%kind /= "split") then
      print *, "No accepted split; check stability and initialization."
      stop
   end if
   if (.not. all(ieee_is_finite(result%x))) error stop "Nonfinite x"
   if (.not. all(ieee_is_finite(result%y))) error stop "Nonfinite y"
   if (.not. ieee_is_finite(result%beta)) error stop "Nonfinite beta"
   if (any(result%x <= 0.0_pr)) error stop "Nonpositive x"
   if (any(result%y <= 0.0_pr)) error stop "Nonpositive y"
   if (result%beta < 0.0_pr .or. result%beta > 1.0_pr) error stop "Invalid beta"

   print *, "beta = ", result%beta
   print *, "x = ", result%x
   print *, "y = ", result%y
   print *, "Composition sums: ", sum(result%x), sum(result%y)
   balance = (1.0_pr-result%beta)*result%x + result%beta*result%y - z
   print *, "Maximum component-balance error: ", maxval(abs(balance))
   call model%lnphi_pt(result%x, P, T, root_type="liquid", lnPhi=lnphi_x)
   call model%lnphi_pt(result%y, P, T, root_type="stable", lnPhi=lnphi_y)
   print *, "Maximum log-fugacity mismatch: ", &
      maxval(abs(log(result%x)+lnphi_x-log(result%y)-lnphi_y))
end program pt_flash
```

The last two diagnostics should be small relative to your required numerical
accuracy. They test balance and phase equilibrium, not global stability. The
parameters here illustrate the API; they are not a fitted mixture data set.

**Current implementation note:** `flash` declares `iters` optional but uses it
as its working iteration counter without a presence check. Pass a declared
integer with `iters=iters` in calls to this routine. Invalid combinations of
specifications can return early without a usable result; validate inputs before
calling rather than relying on `kind` to catch API misuse.

## VT flash

Fragment replacing the call in the PT program above:

```fortran
result = flash(model, z, T, v_spec=1.0_pr, iters=iters)
```

The specified volume is **feed molar volume**, not vessel volume. Divide vessel
volume by total feed moles first. The calculated pressure also satisfies
`v_spec = (1-beta)*Vx + beta*Vy`. As in the PT case, the routine seeks a
two-phase split and is not a general single-phase state classifier.

## Initial K-factors

By default, the flash initializes the K-factors with the Wilson correlation for
`ArModel` models. A better initial estimate can be supplied with `k0 = y/x`:

Fragment (add the import and declaration before executable statements):

```fortran
use yaeos, only: k_wilson

real(pr) :: k0(2)

k0 = k_wilson(model, T=T, P=P)
result = flash(model, z, T, p_spec=P, k0=k0, iters=iters)
```

K-factor initialization is especially important near phase boundaries, for
liquid-liquid splits, and for strongly non-ideal systems.

Note that the flash does not use `k0` unchanged: before the iteration starts it
adjusts the K-factors (scaling them by 1.1 or 0.9) until the Rachford-Rice
equation has a solution with \(0<\beta<1\). This modifies an internal copy,
not the caller's `k0` array. It helps initialize the split equations but does
not establish that the feed is actually two-phase. A genuinely single-phase
feed can be pushed toward a trial split; use stability analysis rather than
treating this adjustment as a phase-classification test.

## `flash_no_beta_limits`

[[flash_no_beta_limits]] is an advanced variant intended for unconstrained
phase fractions. Two caveats apply to the current implementation: the required
`beta0` argument is accepted but never used to initialize the working phase
fraction, so the first Rachford-Rice solve starts from an undefined value; and
the routine prints its iteration history to standard output. Do not rely on
recovery from an uninitialized starting value; the behavior is not defined by
the interface. Signature illustration only (not a
recommended runnable example):

```fortran
result = flash_no_beta_limits( &
   model, z, T, beta0=0.5_pr, p_spec=P, k0=k0, iters=iters &
   )
```

Even in an unconstrained formulation, a fraction outside `[0,1]` is not a
physical phase amount for the specified feed.

# Second-order flash

[[flash_second_order]] calculates a two-phase split at fixed temperature and
pressure by directly minimizing the reduced Gibbs energy with a second-order
Newton method:

```fortran
type(EquilibriumState) function flash_second_order( &
   model, z, T, P_spec, k0, light_root, iters, tol, max_iters &
   ) result(flash)
```

The main arguments are:

| Argument | Meaning |
|:---------|:--------|
| `model` | `ArModel` or `GeModel` |
| `z` | Global mole fractions |
| `T` | Temperature [K] |
| `P_spec` | Pressure [bar] |
| `k0` | Optional initial K-factors |
| `light_root` | `"stable"` (default) or `"liquid"` for `ArModel` |
| `iters` | Iterations used |
| `tol` | Gradient-norm tolerance, default `1e-8` |
| `max_iters` | Maximum iterations, default `200` |

```fortran
use yaeos, only: flash_second_order

result = flash_second_order(model, z, T, P)
```

For an `ArModel`, the heavy phase is always evaluated with the liquid root and
the light phase with the stable root. Use `light_root="liquid"` for a
liquid-liquid split. For a `GeModel`, both phases are described by activity
coefficients.

The second-order method is an alternative when successive substitution stalls,
not a guarantee of convergence or of the globally stable split. Check `kind`,
balances, and stability again. It rejects an excessive final gradient or
identical compositions; failure is reported with `kind="failed"` and `P=-1`.
Use strictly positive normalized feeds: the scaled variables and logarithms
in this implementation divide by each `z(i)`.

# Flash with activity-coefficient models

`flash` accepts [[GeModel]] objects such as [[NRTL]], [[UNIQUAC]], or
[[UNIFAC]]. In that case the pressure and volume cannot be specified, and an
initial guess for the K-factors is mandatory:

Complete liquid-liquid example:

```fortran
program activity_flash
use yaeos, only: pr, UNIFAC, Groups, setup_unifac, EquilibriumState
use yaeos, only: flash, min_tpd
implicit none

integer, parameter :: nc = 3

type(UNIFAC) :: model
type(Groups) :: molecules(nc)
type(EquilibriumState) :: result
real(pr) :: n(nc)
real(pr) :: w(nc)
real(pr) :: T
real(pr) :: P
real(pr) :: mintpd
integer :: iters

! Define the groups of each molecule
molecules(1)%groups_ids = [1, 42]
molecules(1)%number_of_groups = [1, 1]

molecules(2)%groups_ids = [1, 2]
molecules(2)%number_of_groups = [2, 4]

molecules(3)%groups_ids = [16]
molecules(3)%number_of_groups = [1]

model = setup_unifac(molecules)

n = [0.2_pr, 0.7_pr, 0.1_pr]
T = 250.0_pr
P = 1.0_pr

! Find the composition with the minimum tangent-plane distance
call min_tpd(model, n, P=P, T=T, mintpd=mintpd, w=w)

! Use the incipient-phase composition as the initial K-factors
if (mintpd < -1.0e-6_pr) then
   result = flash(model, n, T, k0=w/n, iters=iters)
   if (result%kind == "split") print *, result%x, result%y, result%beta
else
   print *, "No negative TPD found in this search."
end if
end program activity_flash
```

For `GeModel`, `flash` uses the activity coefficients directly and solves only
the isofugacity and mass-balance equations. A stability analysis with
[[min_tpd]] is the recommended initialization strategy. The composition and
phase fraction of the result are meaningful; the pressure and the two phase
volumes are **not** set by the routine and contain undefined values, so do not
read them. Omitting `k0` for a `GeModel` stops the program with an error
message rather than returning a failed result.

# Multiphase flash

[[pt_mp_flash]] is a limited automatic wrapper for `ArModel` systems. Its
`np` counts main phases **excluding** a reference phase; `np=1` represents
two phases, not one. Read [Multiphase flash](multiphase.html) before using it:
the current implementation has incomplete single-phase handling and does not
repeat phase addition until global stability is established.

# Practical notes

- A failed flash means no accepted two-phase result was obtained. It may reflect
  a single-phase feed, poor initialization, a trivial solution, or numerical
  failure. It does **not** establish that the feed is stable.
- Always check `result%kind` before using the compositions.
- Use `k0` from a previous, nearby calculation when building tables or phase
  envelopes. This reduces iterations and avoids convergence to a trivial
  solution.
- For liquid-liquid equilibria, initialize with a stability analysis and
  consider `flash_second_order(..., light_root="liquid")`.
- Never use compositions or volumes from a failed state; fields may be unset
  or contain placeholders. For a `GeModel` the pressure and phase volumes are
  undefined even in a successful result.
- See the runnable examples in `example/basics/3_phase_split.f90` and
  `example/extra/flash.f90`.
