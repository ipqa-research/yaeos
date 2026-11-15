---
title: Cubics
ordered_subpage: alpha.md
                 mixing.md
---

[TOC]

# Cubics

## What a cubic model does

A cubic equation of state is a compact way to estimate fluid pressure, density
and phase behavior from a small set of pure-component data. PR and SRK are
useful starting points for hydrocarbon mixtures; choosing one does not guarantee
accurate liquid densities or strongly polar-mixture behavior. Compare with
measurements in your temperature, pressure and composition range.

Think of the model as a balance between molecules taking up space and molecules
attracting one another. The first term below accounts for excluded volume; the
second reduces pressure because of attraction. Here \(P\) is pressure, \(T\)
is absolute temperature, \(R\) is the gas constant, and \(v\) is **molar** volume.
The covolume \(b\) has units L/mol, the critical attraction parameter \(a_c\)
has units bar L²/mol², and \(\delta_1\), \(\delta_2\) and the temperature
correction \(\alpha\) are dimensionless. Reduced temperature is
\(T_r=T/T_c\), where \(T_c\) is the pure-component critical temperature.
The generic pure-fluid equation [1] is

\[
    P = \frac{RT}{v-b} - \frac{a_c\alpha(T_r)}{(v + \delta_1 b)(v + \delta_2 b)}
\]

For a mixture, the [mixing rule](mixing.html) combines pure-component parameters.
The [alpha function](alpha.html) first adjusts each pure component's attraction
with temperature: these are different jobs. PR, SRK and RKPR default to quadratic
mixing with zero binary corrections unless supplied; PSRK uses an excess-Gibbs
mixing rule instead.

## Required data and a first calculation

For PR and SRK, supply arrays `tc` (K), `pc` (bar), and `w` (dimensionless
acentric factor). The acentric factor summarizes pure-fluid vapor-pressure
behavior; it is not a composition. Every array must use the same component
order. Optional `kij` and `lij` are dimensionless square matrices in that order,
normally symmetric with zero diagonals. They require a suitable parameter source
or fitting; zero is a baseline, not a universal fit.

This complete program uses methane then n-butane. It evaluates a specified
homogeneous state without asking an equilibrium solver to find a state.
The API takes mole **amounts** `n` (mol) and **total** `V` (L), not molar volume.
Here the one-mole basis makes their numerical values coincide.

```fortran
program cubic_pressure
   use yaeos, only: pr, R, CubicEoS, PengRobinson76
   implicit none
   type(CubicEoS) :: model
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: n(2)
   real(pr) :: T
   real(pr) :: V
   real(pr) :: P
   real(pr) :: Pideal

   tc = [190.564_pr, 425.12_pr]
   pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]
   model = PengRobinson76(tc, pc, w)
   n = [0.3_pr, 0.7_pr]
   T = 300.0_pr
   V = 10.0_pr
   call model%pressure(n, V, T, P)
   Pideal = sum(n)*R*T/V
   print *, "Pressure [bar]: ", P
   print *, "Ideal-gas pressure [bar]: ", Pideal
   print *, "Compressibility Z [-]: ", P/Pideal
end program cubic_pressure
```

Link this program with `yaeos` and its dependencies in your configured Fortran/fpm
environment. `pr` is the library's real kind; `_pr` literals avoid rounding inputs
to default real precision. The dimensionless compressibility factor
\(Z=Pv/(RT)\) should approach one at small density. Scaling both `n` and `V`
by the same positive factor should leave pressure unchanged.

## Model variants and their parameters

Below, \(P_c\) is critical pressure and \(\omega\) is acentric factor.
The dimensionless alpha coefficient \(k\) is a pure-component quantity,
**not** a binary interaction \(k_{ij}\). The `b` listed for each model is the
pure-component covolume; for PR and SRK it is also what the quadratic mixing
rule combines.

## SoaveRedlichKwong
Fortran definition: [[SoaveRedlichKwong]]. The SoaveRedlichKwong EoS uses the
following \(\alpha\) function and its respective correlation for \(k\).

\[ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 \]
\[ k = 0.48 + 1.574 \omega - 0.175 \omega^2  \]
\[ a_c = 0.427480  \frac{R^2 T_c^2}{P_c} \]
\[ b = 0.086640  \frac{R T_c}{P_c} \]
\[ \delta_1 = 1 \]
\[ \delta_2 = 0 \]

## PengRobinson76
Fortran definition: [[PengRobinson76]]. The Peng-Robinson EoS uses
the following \(\alpha\) function and correlation for k.

\[ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 \]
\[ k = 0.37464 + 1.54226 \omega - 0.26993 \omega^2 \]
\[ a_c = 0.45723553  \frac{R^2 T_c^2}{P_c} \]
\[ b = 0.07779607 \frac{R T_c}{P_c} \]
\[ \delta_1 = 1 + \sqrt{2} \]
\[ \delta_2 = 1 - \sqrt{2} \]

## PengRobinson78
Fortran definition: [[PengRobinson78]]. The Peng-Robinson (1978) EoS is an
variant of the original Peng-Robinson (1976) equation with a different
high-acentric-factor correlation. It uses the following \(\alpha\) function and
correlation for k.

\[ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 \]

\[
k =
\begin{cases} 
0.37464 + 1.54226 \omega - 0.26992 \omega^2 & \text{if} \quad \omega \leq 0.491 \\
0.379642 + 1.48503 \omega - 0.164423 \omega^2 + 0.016666 \omega^3 & \text{if} \quad \omega > 0.491
\end{cases}
\]

\[ a_c = 0.45723553  \frac{R^2 T_c^2}{P_c} \]

\[ b = 0.07779607  \frac{R T_c}{P_c} \]

\[ \delta_1 = 1 + \sqrt{2} \]

\[ \delta_2 = 1 - \sqrt{2} \]


## RKPR
[[RKPR]]

The RKPR EoS [2] extends the classical formulation of Cubic Equations of State by
freeing the parameter \(\delta_1\) and setting \(\delta_2 =
\frac{1-\delta_1}{1+\delta_1}\) [1]. This extra degree provides extra ways of
implementing the equation in comparison of other Cubic EoS (like PR and SRK)
which are limited to definition of their critical constants.

Besides that extra parameter, the RKPR includes another \(\alpha\)
function:

\[
 \alpha(T_r) = \left(\frac{3}{2+T_r}\right)^k
\]

Its constructor is `RKPR(tc, pc, w, zc, kij, lij, delta_1, k)`; arguments after
`zc` are optional. Here `zc` is the experimental critical compressibility
\(Z_c=P_c v_c/(RT_c)\), with \(v_c\) the critical molar volume. RKPR needs
this additional pure-component information. Optional `delta_1` and `k` arrays
each have one entry per component.

The implementation uses an effective critical compressibility `1.168*zc` in
its default shape correlation and numerically readjusts unspecified `k` values
to match the acentric-factor definition at \(T_r=0.7\). This numerical fit is
not an unconditional convergence guarantee. The implemented correlations are
documented with [[RKPR]]; do not substitute experimental `zc` directly into an
unrelated published correlation. Check that correlated parameters are finite,
particularly near the fractional-power correlation's domain boundary.

Because the \(k\) refit solves a saturation condition per component, building
an RKPR model is noticeably more expensive than building PR or SRK, and it can
be sensitive to the input `zc` values.


## PSRK
[[PSRK]]

The PSRK EoS [3,4] (Predictive Soave-Redlich-Kwong) is a Soave-Redlich-Kwong
cubic equation of state that uses the [[MHV]] mixrule and a PSRK-UNIFAC model.
The alpha function is Soave by default, or [[AlphaMathiasCopeman]] when the
optional `c1`, `c2`, and `c3` parameters are provided.

PSRK additionally needs one group decomposition per component, using the
PSRK-UNIFAC group definitions and available interaction parameters. It does not
infer molecular structure from critical constants. Supply all three alpha
coefficient arrays together if overriding its default. See the
[group-contribution guide](../../excessmodels/unifac_models/index.html).

## Interpreting results and avoiding false confidence

- Use positive temperatures and volumes, nonnegative amounts and a positive
  total amount. Volume must exceed the mixture covolume; approaching that
  boundary makes pressure diverge.
- `pressure` evaluates a mathematical homogeneous state. It does not decide
  whether that state will separate into liquid and vapor.
- At fixed pressure and temperature there can be several volume roots.
  Requesting a liquid or vapor root is a branch choice, not proof that this
  phase is globally stable. Negative pressures on stretched branches are not
  a predicted operating pressure for an equilibrium vessel.
- A converged flash cannot repair unsuitable parameters. Conversely, a failed
  flash is not proof the EOS is invalid. Check units, root selection, initial
  guesses and material balances separately.

## Next reading

Continue with [mixing rules](mixing.html), [alpha functions](alpha.html), and
[phase equilibrium](../../phase_equilibrium/index.html). The
[EOS property guide](../index.html) explains the shared API, and
[thermodynamic fundamentals](../../fundamentals.html) covers amounts,
compositions, units and the difference between a model and a solver.


# References
[1] Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
Fundamentals & computational aspects (2. ed). Tie-Line Publications.

[2] M. Cismondi, J. Mollerup, Development and application of a three-parameter 
RK–PR equation of state, Fluid Phase Equilibria 232 (2005) 74–89.
https://doi.org/10.1016/j.fluid.2005.03.020.

[3] Holderbaum, T., & Gmehling, J. (1991). PSRK: A Group Contribution Equation
of State Based on UNIFAC. Fluid Phase Equilibria, 70(2-3), 251-265.
https://doi.org/10.1016/0378-3812(91)85038-V

[4] Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., & Gmehling, J.
(2005). PSRK group contribution equation of state: Comprehensive revision and
extension IV, including critical constants and α-function parameters for 1000
components. Fluid Phase Equilibria, 227(2), 157-164.
https://doi.org/10.1016/j.fluid.2004.11.002
