---
title: Critical points
---

[TOC]

# Critical points

A critical point is the state where two coexisting phases become
indistinguishable: their compositions, molar volumes, and densities converge,
and the boundary between them disappears. For a pure substance it is a single
point in pressure-temperature space. For a mixture it also depends on
composition, so the critical states form a **critical line** in
pressure-temperature-composition space.

Near a critical point the usual "which phase is this?" questions stop having
sharp answers. That is why critical calculations are useful: they mark the end
of a phase envelope, and they explain why a flash or a saturation solver becomes
ill-conditioned in that region.

`yaeos` provides:

- [[critical_point]]: a single critical point at a specified composition.
- [[critical_line]]: a critical line between two limiting compositions.
- [[get_critical_constants]]: pure-component critical constants of a model.

# Theory

For a pure fluid, the critical point is the point where the pressure-volume
curve develops an inflection:

$$
   \left(\frac{\partial P}{\partial V}\right)_T = 0, \qquad
   \left(\frac{\partial^2 P}{\partial V^2}\right)_T = 0
$$

For a mixture, composition fluctuations matter as well as changes in density;
the pure-fluid pressure-inflection test is not a general mixture criticality
criterion. Write the total Helmholtz energy as an ideal-gas contribution plus
a residual (real-fluid) contribution,

$$
   A(\mathbf n,V,T) = A^{\mathrm{ideal}}(\mathbf n,V,T) + A^{r}(\mathbf n,V,T).
$$

The mixture loses stability when the matrix of second derivatives of the
**total** Helmholtz energy with respect to the component mole numbers becomes
singular:

$$
   \left(\frac{\partial^2 A}{\partial n_i \partial n_j}\right)_{T,V}
   \text{ is singular}.
$$

The ideal contribution is not optional. In the Hessian divided by \(RT\), it
supplies the positive diagonal elements \(1/n_i\). An ideal gas with positive
component amounts is stable against these fluctuations: this term does not
make it unstable at low temperature. In a real fluid, the residual contribution
can offset that stabilizing curvature. Using the residual Hessian alone would
therefore test a different, generally incorrect condition.

This is how the implementation proceeds. The fugacity derivative returned by
`lnfug_vt` satisfies

$$
   \frac{\partial \ln f_i}{\partial n_j}
   = \frac{1}{n_i}\delta_{ij} + \frac{1}{RT}\frac{\partial^2 A^{r}}{\partial n_i \partial n_j},
$$

that is, the ideal-mixing term is already included. The stability matrix built
from it is

$$
   M_{ij} = \sqrt{z_i z_j}\;\frac{\partial \ln f_i}{\partial n_j},
$$

and the implementation selects the eigenvalue closest to zero in absolute
value, denoted \(\lambda_1\) here. This selection is an implementation choice;
it is not by itself a test that all other eigenvalues are nonnegative.
The implementation evaluates \(\lambda_1\) at
the state and along a small composition perturbation
\(\mathbf n = \mathbf z + s\,\mathbf u\sqrt{\mathbf z}\), where \(\mathbf u\) is
the eigenvector associated with the selected eigenvalue. The solver imposes

$$
   \lambda_1(0) = 0
   \qquad\text{and}\qquad
   \frac{d\lambda_1}{ds}\Big|_{s=0} = 0,
$$

plus the pressure identity \(P(\mathbf n,V,T) = \exp(X_4)\) and the requested
specification. The second condition is evaluated by a central difference with
\(\varepsilon = 10^{-4}\), and the whole system is solved by Newton iteration.

# Critical-point specification

Critical-point calculations solve a system of equations involving the
stability limit and the criticality conditions. One variable must be
specified. The `spec_CP` object provides the allowed specifications:

| Field | Meaning |
|:------|:--------|
| `spec_CP%a` | Composition parameter \(\alpha\) |
| `spec_CP%V` | Volume |
| `spec_CP%T` | Temperature |
| `spec_CP%P` | Pressure |

The composition of the mixture is defined as:

$$
    z = \alpha z_i + (1 - \alpha) z_0
$$

where \(z_0\) and \(z_i\) are the two limiting compositions. With
\(z_0=[1,0]\) and \(z_i=[0,1]\), \(\alpha\) is simply the mole fraction of the
second component. With more general limiting compositions, \(\alpha\) is a
mixing parameter and the actual mole fractions are the linear combination above.

# `critical_point`

```fortran
type(EquilibriumState) function critical_point( &
   model, z0, zi, spec, S, max_iters, V0, T0, a0, P0 &
   )
   class(ArModel), intent(in) :: model
   real(pr), intent(in) :: z0(:)   !! First limiting composition
   real(pr), intent(in) :: zi(:)   !! Second limiting composition
   integer, intent(in) :: spec      !! spec_CP field
   real(pr), intent(in) :: S        !! Specification value
   integer, intent(in) :: max_iters
   real(pr), optional, intent(in) :: V0
   real(pr), optional, intent(in) :: T0
   real(pr), optional, intent(in) :: a0
   real(pr), optional, intent(in) :: P0
end function critical_point
```

The specification value `S` follows the natural variables of the solver,
\((\alpha, \ln V, \ln T, \ln P)\):

- For `spec_CP%a`, `S` is the composition parameter \(\alpha\) itself, in
  \([0,1]\).
- For `spec_CP%T`, `S` is \(\ln T\), with \(T\) in kelvin.
- For `spec_CP%P`, `S` is \(\ln P\), with \(P\) in bar.
- For `spec_CP%V`, `S` is \(\ln V\), with \(V\) in L/mol.

The `V0`, `T0`, `P0`, and `a0` arguments are initial guesses, given in the
ordinary units (K, bar, L/mol, and \(\alpha\) respectively), not as logarithms.
Only the specified variable `S` uses the logarithmic convention. Passing a
temperature in kelvin as `S` for `spec_CP%T` therefore produces a nonsense
state rather than an error.

## Example

```fortran
program critical_example
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState
   use yaeos, only: critical_point, spec_CP
   implicit none

   type(CubicEoS) :: model
   type(EquilibriumState) :: cp
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: z0(2)
   real(pr) :: zi(2)
   real(pr) :: a

   ! Carbon dioxide + n-hexane
   tc = [304.21_pr, 507.6_pr]
   pc = [73.83_pr, 30.25_pr]
   w = [0.223621_pr, 0.301261_pr]

   model = PengRobinson76(tc, pc, w)

   z0 = [1.0_pr, 0.0_pr]  ! Pure carbon dioxide
   zi = [0.0_pr, 1.0_pr]  ! Pure n-hexane
   a = 0.5_pr             ! Equimolar mixture

   cp = critical_point( &
      model, z0, zi, spec=spec_CP%a, S=a, max_iters=1000, a0=a &
      )

   print *, "Critical temperature: ", cp%T, " K"
   print *, "Critical pressure: ", cp%P, " bar"
   print *, "Critical volume: ", cp%Vx, " L/mol"
   print *, "Critical composition: ", cp%x
   print *, "Iterations: ", cp%iters
end program critical_example
```

The result is an [[EquilibriumState]] with `kind="critical"`, `x == y == z`,
and `Vx == Vy`. `beta` is zero. With the parameters above the calculation
converges in a few tens of Newton iterations to about 461.2 K and 76.4 bar.

`critical_point` has no separate failure flag: `kind` is set to `"critical"`
when the iteration ends, whether or not the residuals are small. Check the
returned state instead: the pressure must be positive and finite, the
composition must be normalized and nonnegative, and the result must be
consistent with an independent calculation.

# `critical_line`

`critical_line` continues a critical point along the composition parameter
\(\alpha\). It uses the same continuation method as the phase envelopes.

```fortran
type(CriticalLine) function critical_line( &
   model, a0, z0, zi, ns0, S0, dS0, &
   v0, t0, p0, max_points, maxP, first_point, stability_analysis &
   )
   class(ArModel), intent(in) :: model
   real(pr), intent(in) :: a0       !! Initial alpha
   real(pr), intent(in) :: z0(:)
   real(pr), intent(in) :: zi(:)
   integer, intent(in) :: ns0       !! Specification variable
   real(pr), intent(in) :: S0       !! Initial specification value
   real(pr), intent(in) :: dS0      !! Initial step
   real(pr), optional, intent(in) :: v0
   real(pr), optional, intent(in) :: t0
   real(pr), optional, intent(in) :: p0
   integer, optional, intent(in) :: max_points
   real(pr), optional, intent(in) :: maxP
   type(EquilibriumState), optional, intent(in) :: first_point
   logical, optional :: stability_analysis
end function critical_line
```

`ns0` is usually `spec_CP%a`. In that case `S0` is the initial \(\alpha\) and
`dS0` is the step in \(\alpha\). A negative `dS0` scans in the opposite
direction. The optional `first_point` allows the continuation to start from an
already converged critical point; without it the routine starts from an
initial guess built from the component critical constants and a volume
estimate.

When `ns0` is not `spec_CP%a`, `S0` follows the same logarithmic convention as
`critical_point`: `log(T)`, `log(P)`, or `log(V)`.

The optional `first_point` is used for two things: the volume, temperature and
pressure guesses, and the composition parameter, which is taken from
`first_point%x(2)`. That second element is the mole fraction of the second
component only for a binary mixture whose limiting compositions are the pure
components. For other mixtures, set the starting point with `a0` and `S0`
instead.

```fortran
use yaeos, only: CriticalLine, critical_line, critical_point, spec_CP

type(CriticalLine) :: cl
type(EquilibriumState) :: cp
real(pr) :: a

a = 0.5_pr

cp = critical_point(model, z0, zi, spec=spec_CP%a, S=a, max_iters=1000, a0=a)

cl = critical_line( &
   model, a0=a, z0=z0, zi=zi, &
   ns0=spec_CP%a, S0=a, dS0=0.01_pr, &
   first_point=cp, max_points=1000 &
   )

print *, "Number of points: ", size(cl%a)
print *, "Alpha range: ", cl%a(1), cl%a(size(cl%a))
```

The `CriticalLine` type contains:

| Field | Meaning |
|:------|:--------|
| `a` | Composition parameter \(\alpha\) at each point |
| `z0`, `zi` | Limiting compositions |
| `P` | Pressure [bar] |
| `V` | Molar volume [L/mol] |
| `T` | Temperature [K] |
| `ns` | Specification variable used at each point |
| `iters` | Iterations used at each point |
| `CEP` | Critical end point, when one is found |

The number of points is decided by the continuation itself: it stops when the
line turns, when a point fails to converge, when the pressure exceeds `maxP`
(default 2500 bar), or when `max_points` is reached. An empty line is a
legitimate outcome and means the starting point was not usable, not that the
mixture has no critical behaviour. Starting from a converged critical point
matters: with the example above, a continuation started from a poor point
returned an empty line.

`CEP` is only filled when `stability_analysis` is requested and a critical end
point is actually detected. The detection routine is limited to binary
mixtures and scans a fixed grid of trial compositions, so a null `CEP` is not
evidence that no critical end point exists.

# `get_critical_constants`

[[get_critical_constants]] calculates the pure-component critical constants of
an `ArModel` and updates the model's `components` object in place:

```fortran
use yaeos__equilibria_critical, only: get_critical_constants

call get_critical_constants(model)

! The critical constants are now available in the model
print *, model%components%Tc
print *, model%components%Pc
print *, model%components%w
```

This routine is **not** re-exported by the top-level `yaeos` module; it must be
imported from the implementation module `yaeos__equilibria_critical` as shown
above. Use it only with a model that has already allocated `Tc`, `Pc`, and `w`
for its components.

The routine solves a critical point for each pure component, then estimates
the acentric factor from the saturation pressure at \(0.7\,T_c\):

$$
   \omega = -1 - \log_{10}\left(\frac{P^{\mathrm{sat}}(0.7\,T_c)}{P_c}\right).
$$

It is used internally by models such as PC-SAFT, where the critical constants
are not known analytically and must be obtained from the equation of state
itself.

# Practical notes

- Critical-point calculations are sensitive to the initial guesses. Use `a0`,
  `T0`, `P0`, and `V0` when they are available.
- Specify temperature or pressure in preference to volume. A run with
  `spec_CP%V` for the example mixture returned a state with a pressure of
  order \(10^{-7}\) bar, which is a converged iteration of the wrong problem.
- The specification value is logarithmic for `spec_CP%T`, `spec_CP%P`, and
  `spec_CP%V`; the initial guesses are not.
- `critical_point` does not report convergence failure. Validate the returned
  temperature, pressure, and composition, and compare with a second
  specification when possible.
- A critical line can cross a critical end point or a liquid-liquid region.
  Check the `CEP` field and the stability of the calculated points.
- The critical line and the phase envelope should meet at the critical point.
  Use this as a consistency check.
- Near a critical point the two phases are nearly identical, so the
  saturation and flash routines lose their conditioning: phase fractions,
  volumes, and compositions all become uncertain at the same time. Treat
  results close to a critical point as qualitative.
- Critical calculations are expensive because they require second and third
  derivatives of the Helmholtz energy through numerical differentiation of the
  criticality conditions.
