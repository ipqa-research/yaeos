---
title: Alpha functions
---

[TOC]

# Why a cubic EOS needs an alpha function

Molecular attractions do not have the same effect at every temperature. A cubic
EOS represents that change by multiplying a reference attraction parameter by
a dimensionless temperature correction. This is the **pure-component** part of
the model: a [mixing rule](mixing.html) subsequently combines components. Changing
alpha is not a substitute for fitting binary interactions.

For component \(i\), let \(T\) be absolute temperature in K, \(T_{c,i}\) its
critical temperature in K, and \(T_{r,i}=T/T_{c,i}\) reduced temperature.
The critical attraction parameter \(a_{c,i}\) and the temperature-dependent
parameter \(a_i(T)\) both have units bar L²/mol². Then

\[
a_i(T)=a_{c,i}\alpha_i(T_{r,i}).
\]

Built-in cubic constructors already provide an alpha function. Keep that default
unless you have pure-fluid vapor-pressure data or a compatible published fit
that justifies a change. A more flexible alpha can improve a fit over its
calibration interval without ensuring good extrapolation or mixture predictions.

## Implemented choices

Below, \(k\), \(c_1\), \(c_2\) and \(c_3\) are dimensionless **pure-component**
coefficients. They are not the binary interaction \(k_{ij}\), nor the NRTL
non-randomness parameter also often called alpha. Each coefficient array has
one entry per component in exactly the constructor's component order.

### Soave

[[AlphaSoave]] uses

\[
\alpha(T_r)=[1+k(1-\sqrt{T_r})]^2.
\]

With \(\omega\) the dimensionless acentric factor, the constructors use:

| Model | \(k(\omega)\) |
|:------|:--------------|
| Peng-Robinson 1976 | \(0.37464+1.54226\omega-0.26993\omega^2\) |
| Peng-Robinson 1978, \(\omega\leq0.491\) | \(0.37464+1.54226\omega-0.26992\omega^2\) |
| Peng-Robinson 1978, \(\omega>0.491\) | \(0.379642+1.48503\omega-0.164423\omega^2+0.016666\omega^3\) |
| Soave-Redlich-Kwong | \(0.48+1.574\omega-0.175\omega^2\) |

Its reduced-temperature derivatives are

\[
\frac{d\alpha}{dT_r}=\frac{k[k(\sqrt{T_r}-1)-1]}{\sqrt{T_r}},\qquad
\frac{d^2\alpha}{dT_r^2}=\frac{k(k+1)}{2T_r^{3/2}}.
\]

### RKPR

[[AlphaRKPR]] uses its component-specific exponent \(k\):

\[
\alpha(T_r)=\left(\frac{3}{2+T_r}\right)^k,\qquad
\frac{d\alpha}{dT_r}=-\frac{k\alpha}{2+T_r},\qquad
\frac{d^2\alpha}{dT_r^2}=-\frac{k+1}{2+T_r}\frac{d\alpha}{dT_r}.
\]

The [[RKPR]] constructor estimates and numerically refits unspecified exponents
using the acentric factor and critical-compressibility information. See the
[cubic overview](index.html) for its additional inputs and initialization caveats.

### Mathias-Copeman

[[AlphaMathiasCopeman]] adds flexibility below the critical temperature.
Define \(u=1-\sqrt{T_r}\). Then

\[
\alpha(T_r)=
\begin{cases}
(1+c_1u+c_2u^2+c_3u^3)^2,&T_r\leq1,\\
(1+c_1u)^2,&T_r>1.
\end{cases}
\]

Thus `c2` and `c3` affect only the subcritical branch. [[PSRK]] uses this
function when all three optional coefficient arrays are supplied; otherwise
it uses Soave. The implementation supplies analytical derivatives in
`src/models/residual_helmholtz/cubic/alphas/alphas.f90`.

## Complete example: install and inspect an alpha

This program installs the built-in Mathias-Copeman type rather than reimplementing
it. The methane/n-butane critical data have the same order as the coefficient
arrays. **The coefficients are illustrative:** setting `c2=c3=0` and using PR76's
`k` as `c1` reproduces the Soave form, providing a useful consistency check,
not a new fitted parameter set.

```fortran
program alpha_comparison
   use yaeos, only: pr, CubicEoS, PengRobinson76, AlphaMathiasCopeman
   implicit none
   type(CubicEoS) :: model
   type(AlphaMathiasCopeman) :: correction
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: n(2)
   real(pr) :: Tr(2)
   real(pr) :: value(2)
   real(pr) :: first(2)
   real(pr) :: second(2)
   real(pr) :: T
   real(pr) :: V
   real(pr) :: Pbefore
   real(pr) :: Pafter

   tc = [190.564_pr, 425.12_pr]
   pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]
   model = PengRobinson76(tc, pc, w)
   n = [0.3_pr, 0.7_pr] ! mol
   T = 300.0_pr ! K
   V = 10.0_pr ! total L
   call model%pressure(n, V, T, Pbefore)

   correction%c1 = 0.37464_pr + 1.54226_pr*w - 0.26993_pr*w**2
   correction%c2 = [0.0_pr, 0.0_pr]
   correction%c3 = [0.0_pr, 0.0_pr]
   call model%set_alpha(correction)
   call model%pressure(n, V, T, Pafter)
   Tr = T/tc
   call correction%alpha(Tr, value, first, second)
   print *, "Alpha values [-]: ", value
   print *, "d(alpha)/d(Tr): ", first
   print *, "d2(alpha)/d(Tr)2: ", second
   print *, "Pressure difference [bar]: ", Pafter-Pbefore
end program alpha_comparison
```

The pressure difference should be zero to floating-point roundoff. Check each
component's alpha at `Tr=1`: these built-in forms give one there. At 300 K the
two components in this example are on different sides of their critical
temperatures, so both Mathias-Copeman branches are exercised.

`set_alpha` deallocates the previous polymorphic component before assigning the
new one. It is the preferred replacement method. To investigate a real fit,
replace the illustrative coefficients with documented values fitted to this
EOS, then compare calculated pure-fluid saturation pressures with data.

## Custom functions and derivative checks

[[AlphaFunction]] is the abstract parent. A new type extends it and overrides
`alpha`. This **interface sketch**, not a standalone program, shows the
required arguments:

```fortran
subroutine alpha(self, Tr, a, dadt, dadt2)
   class(MyAlpha), intent(in) :: self
   real(pr), intent(in) :: Tr(:)
   real(pr), intent(out) :: a(:)
   real(pr), intent(out) :: dadt(:)
   real(pr), intent(out) :: dadt2(:)
   ! Supply alpha and its first two derivatives with respect to Tr.
end subroutine alpha
```

All arrays have one entry per component. Despite the argument name `a`, this
procedure returns dimensionless alpha, **not** the attraction parameter.
Despite the names `dadt` and `dadt2`, these are derivatives with respect to
**reduced** temperature. The cubic model performs the chain-rule conversion:

\[
\frac{da_i}{dT}=\frac{a_{c,i}}{T_{c,i}}\frac{d\alpha_i}{dT_{r,i}},\qquad
\frac{d^2a_i}{dT^2}=\frac{a_{c,i}}{T_{c,i}^2}
\frac{d^2\alpha_i}{dT_{r,i}^2}.
\]

Returning absolute-temperature derivatives here would apply the conversion
twice. Pressure alone may not reveal the error; residual energy properties and
solver derivatives can be wrong. Compare analytical derivatives with finite
differences away from branch boundaries, and inspect continuity at boundaries
separately. Use positive `Tr`; square roots and inverse powers are undefined
or singular at nonpositive values.

Squared expressions ensure nonnegative alpha, but do not guarantee physically
reasonable high-temperature trends. Check the fitted range, critical-point
normalization and derivative behavior before extrapolating. Numerical
convergence of a saturation calculation is distinct from accuracy of the fit.

## Next reading

`example/advanced/new_alpha_function.f90` contains a complete custom-type
implementation. Continue with [mixing rules](mixing.html), the
[EOS property guide](../index.html),
[thermodynamic fundamentals](../../fundamentals.html), and
[saturation calculations](../../phase_equilibrium/saturation_points.html).
Model references are collected on the [cubic overview](index.html).
