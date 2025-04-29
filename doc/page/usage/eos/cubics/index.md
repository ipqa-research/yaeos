---
title: Cubics
---

[TOC]

All our Cubic Equations of State are implemented based on the generic Cubic
Equation [1]:

$$
    P = \frac{RT}{V-b} - \frac{a_c\alpha(T_r)}{(V + \delta_1 b)(V - \delta_2 b)}
$$

All cubic equation of state allows to include the \(k_{ij}\) and \(l_{ij}\)
matrices. Using by default Classic Van der Waals mixing rules. For more
information about mixing rules look at [Mixing Rules](mixing.html). Moreover,
cubic equations of state allows to modify their default \(\alpha\) function,
look at [Alpha functions](alpha.html).


## SoaveRedlichKwong
[[SoaveRedlichKwong]]

$$ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 $$
$$ k = 0.48 + 1.574 \omega - 0.175 \omega^2  $$
$$ a_c = 0.427480  \frac{R^2 T_c^2}{P_c} $$
$$ b = 0.086640  \frac{R T_c}{P_c} $$
$$ \delta_1 = 1 $$
$$ \delta_2 = 0 $$

### Example

```fortran
use yaeos, only: pr, SoaveRedlichKwong, ArModel

class(ArModel), allocatable :: model

real(pr) :: tc(2), pc(2), w(2)
real(pr) :: n(2), p, t, v

Tc = [190.564, 425.12]              ! Critical temperatures [K]
Pc = [45.99, 37.96]                 ! Critical pressures [bar]
w = [0.0115478, 0.200164]           ! Acentric factors

model = SoaveRedlichKwong(Tc, Pc, w)

n = [4.0, 6.0]                      ! Composition [mol]
p = 10.0                            ! Pressure [bar]
t = 300.0                           ! Temperature [K]

call model%volume(n, p, t, v, root_type="stable")

print *, "Volume: ", v, "L"
```

## PengRobinson76
[[PengRobinson76]]

$$ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 $$
$$ k = 0.37464 + 1.54226 \omega - 0.26993 \omega^2 $$
$$ a_c = 0.45723553  \frac{R^2 T_c^2}{P_c} $$
$$ b = 0.07779607 \frac{R T_c}{P_c} $$
$$ \delta_1 = 1 + \sqrt{2} $$
$$ \delta_2 = 1 - \sqrt{2} $$

### Example

```fortran
use yaeos, only: pr, PengRobinson76, ArModel

class(ArModel), allocatable :: model

real(pr) :: tc(2), pc(2), w(2)
real(pr) :: n(2), p, t, v

Tc = [190.564, 425.12]              ! Critical temperatures [K]
Pc = [45.99, 37.96]                 ! Critical pressures [bar]
w = [0.0115478, 0.200164]           ! Acentric factors

model = PengRobinson76(Tc, Pc, w)

n = [4.0, 6.0]                      ! Composition [mol]
p = 10.0                            ! Pressure [bar]
t = 300.0                           ! Temperature [K]

call model%volume(n, p, t, v, root_type="stable")

print *, "Volume: ", v, "L"
```

## PengRobinson78
[[PengRobinson78]]

$$ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 $$
$$
k =
\begin{cases} 
0.37464 + 1.54226 \omega - 0.26992 \omega^2 & \text{if} \quad \omega \leq 0.491 \\
0.37464 + 1.48503 \omega - 0.16442 \omega^2 + 0.016666 \omega^3 & \text{if} \quad \omega > 0.491
\end{cases}
$$
$$ a_c = 0.45723553  \frac{R^2 T_c^2}{P_c} $$
$$ b = 0.07779607  \frac{R T_c}{P_c} $$
$$ \delta_1 = 1 + \sqrt{2} $$
$$ \delta_2 = 1 - \sqrt{2} $$

### Example

```fortran
use yaeos, only: pr, PengRobinson78, ArModel

class(ArModel), allocatable :: model

real(pr) :: tc(2), pc(2), w(2)
real(pr) :: n(2), p, t, v

Tc = [190.564, 425.12]              ! Critical temperatures [K]
Pc = [45.99, 37.96]                 ! Critical pressures [bar]
w = [0.0115478, 0.200164]           ! Acentric factors

model = PengRobinson78(Tc, Pc, w)

n = [4.0, 6.0]                      ! Composition [mol]
p = 10.0                            ! Pressure [bar]
t = 300.0                           ! Temperature [K]

call model%volume(n, p, t, v, root_type="stable")

print *, "Volume: ", v, "L"
```


## RKPR
[[RKPR]]

The RKPR EoS extends the classical formulation of Cubic Equations of State by
freeing the parameter \(\delta_1\) and setting \(\delta_2 =
\frac{1+\delta_1}{1-\delta_1}\) [1]. This extra degree provides extra ways of
implementing the equation in comparison of other Cubic EoS (like PR and SRK)
which are limited to definition of their critical constants.

Besides that extra parameter, the RKRR includes another \(\alpha\)
function:

$$
 \alpha(T_r) = \left(\frac{3}{2+T_r}\right)^k
$$

In this implementation we take the simplest form which correlates
the extra parameter to the critical compressibility factor \(Z_c\) and
the \(k\) parameter of the \(\alpha\) function to \(Z_c\) and \(\omega\):

$$ \delta_1 = d_1 + d_2 (d_3 - Z_c)^d_4 + d_5 (d_3 - Z_c) ^ d_6 $$
$$ k = (A_1  Z_c + A_0)\omega^2 + (B_1 Z_c + B_0)\omega + (C_1 Z_c + C_0) $$

It is also possible to include the parameters as optional arguments.

### Example

```fortran
use yaeos, only: pr, RKPR, ArModel

class(ArModel), allocatable :: model

real(pr) :: Tc(2), Pc(2), Zc(2), w(2)
real(pr) :: n(2), p, t, v

Tc = [190.564, 425.12]              ! Critical temperatures [K]
Pc = [45.99, 37.96]                 ! Critical pressures [bar]
Zc = [0.286, 0.274]                 ! Critical compressibility factors
w = [0.0115478, 0.200164]           ! Acentric factors

model = RKPR(Tc, Pc, w, Zc)

n = [4.0, 6.0]                      ! Composition [mol]
p = 10.0                            ! Pressure [bar]
t = 300.0                           ! Temperature [K]

call model%volume(n, p, t, v, root_type="stable")

print *, "Volume: ", v, "L"
```


## PSRK
[[PSRK]]

The PSRK EoS [2,3] (Predictive Soave-Redlich-Kwong) is a Soave-Redlich-Kwong
cubic equation of state that uses the [[MHV]] mixrule and the
[[AlphaMathiasCopeman]] \(\alpha\) function. With the MHV mixrule, a Gibbs
excess energy model is required, in the case of PSRK, the PSRK-UNIFAC model is
used.

In the next example of a Methane-Butane mixture, the Mathias-Copeman 
coefficients are given. If the Mathias-Copeman coefficients are not given, the
[[AlphaSoave]] function is used with the provided acentric factors.

The functional groups must be specified from the subgroups ids and its 
occurences. To check the full list of funcitonal groups of the PSRK-UNIFAC 
check the [PSRK-UNIFAC](../../excessmodels/psrk-unifac.html) page.

### Example

```fortran
use yaeos, only: pr, PSRK, ArModel, Groups

class(ArModel), allocatable :: model
type(Groups) :: molecules(2)

real(pr) :: Tc(2), Pc(2), w(2), c1(2), c2(2), c3(2)
real(pr) :: n(2), p, t, v

Tc = [190.564, 425.12]              ! Critical temperatures [K]
Pc = [45.99, 37.96]                 ! Critical pressures [bar]
w = [0.0115478, 0.200164]           ! Acentric factors
c1 = [0.49258, 0.84209]             ! Mathias-Copeman coefficient c1
c2 = [0.0, -0.46406]                ! Mathias-Copeman coefficient c2
c3 = [0.0, 0.84619]                 ! Mathias-Copeman coefficient c3

molecules(1)%groups_ids = [118]
molecules(1)%number_of_groups = [1]

molecules(2)%groups_ids = [1, 2]
molecules(2)%number_of_groups = [2, 2]

model = PSRK(Tc, Pc, w, molecules, c1, c2, c3)

n = [4.0, 6.0]                      ! Composition [mol]
p = 10.0                            ! Pressure [bar]
t = 300.0                           ! Temperature [K]

call model%volume(n, p, t, v, root_type="stable")

print *, "Volume: ", v, "L"
```

## References
[1] Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
Fundamentals & computational aspects (2. ed). Tie-Line Publications.

[2] Holderbaum, T., & Gmehling, J. (1991). PSRK: A Group Contribution Equation
of State Based on UNIFAC. Fluid Phase Equilibria, 70(2-3), 251-265.
https://doi.org/10.1016/0378-3812(91)85038-V

[3] Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., & Gmehling, J.
(2005). PSRK group contribution equation of state: Comprehensive revision and
extension IV, including critical constants and α-function parameters for 1000
components. Fluid Phase Equilibria, 227(2), 157-164.
https://doi.org/10.1016/j.fluid.2004.11.002


