---
title: Cubics
---

[TOC]

Our Cubic Equations of State are implemented based on the generic Cubic
Equation [1]:

\[
    P = \frac{RT}{V-b} - \frac{a_c\alpha(T_r)}{(V + \delta_1 b)(V - \delta_2 b)}
\]

Cubic equation of state allows to include the \(k_{ij}\) and \(l_{ij}\)
matrices. 
Using by default Classic Van der Waals mixing rules.
For more information about mixing rules look at [Mixing Rules](mixing.md).
Moreover, cubic equations of state allows to modify their default 
\(\alpha\) function, look at [Alpha functions](alpha.md).


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
Fortran definition: [[PengRobinson78]]. The Peng-Robinson 78 EoS is an
improved version of the original PengRobinson equation for heavier components. 
This equation ensures a monotonically increasing \(k\) as the values of 
\(\omega\) increasises.
It uses the following \(\alpha\) function and correlation for k.

\[ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 \]
\[
k =
\begin{cases} 
0.37464 + 1.54226 \omega - 0.26992 \omega^2 & \text{if} \quad \omega \leq 0.491 \\
0.37464 + 1.48503 \omega - 0.16442 \omega^2 + 0.016666 \omega^3 & \text{if} \quad \omega > 0.491
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
\frac{1+\delta_1}{1-\delta_1}\) [1]. This extra degree provides extra ways of
implementing the equation in comparison of other Cubic EoS (like PR and SRK)
which are limited to definition of their critical constants.

Besides that extra parameter, the RKRR includes another \(\alpha\)
function:

\[
 \alpha(T_r) = \left(\frac{3}{2+T_r}\right)^k
\]

These two extra parameters can be provided as arguments. But, if they are
not provided they will be calculated by the following correlations:

\[ \delta_1 = d_1 + d_2 (d_3 - Z_c)^d_4 + d_5 (d_3 - Z_c) ^ d_6 \]
\[ k = (A_1  Z_c + A_0)\omega^2 + (B_1 Z_c + B_0)\omega + (C_1 Z_c + C_0) \]

In this implementation, the \(k\) constants (if not provided) will also 
be readjusted to assure that \(\omega = -log_{10}\left(P_r^{sat}\right) - 1\) 
at \(T_r = 0.7\).

It is also possible to include the parameters as optional arguments.


## PSRK
[[PSRK]]

The PSRK EoS [3,4] (Predictive Soave-Redlich-Kwong) is a Soave-Redlich-Kwong
cubic equation of state that uses the [[MHV]] mixrule and the
[[AlphaMathiasCopeman]] \(\alpha\) function. With the MHV mixrule, a Gibbs
excess energy model is required, in the case of PSRK, the PSRK-UNIFAC model is
used.


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
