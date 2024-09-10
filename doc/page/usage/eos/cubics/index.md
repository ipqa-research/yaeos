---
title: Cubics
---

All our Cubic Equations of State are implemented based on the generic Cubic
Equation:

\[
    P = \frac{RT}{V-b} - \frac{a_c\alpha(T_r)}{(V + \delta_1 b)(V - \delta_2 b)}
\]

## SoaveRedlichKwong
[[SoaveRedlichKwong]]

Using the critical constants setup the parameters to use the 
SoaveRedlichKwong Equation of State

- \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
- \[k = 0.48 + 1.574 \omega - 0.175 \omega^2 \]
- \[a_c = 0.427480  R^2 * T_c^2/P_c\]
- \[b = 0.086640  R T_c/P_c\]
- \[\delta_1 = 1\]
- \[\delta_2 = 0\]

There is also the optional posibility to include the k_{ij} and l_{ij}
matrices. Using by default Classic Van der Waals mixing rules. For more information
about mixing rules look at [Mixing Rules](mixing.md)

## PengRobinson76
[[PengRobinson76]]

- \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
- \[k = 0.37464 + 1.54226 * \omega - 0.26993 \omega^2 \]
- \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
- \[b = 0.07779607r  R T_c/P_c\]
- \[\delta_1 = 1 + \sqrt{2}\]
- \[\delta_2 = 1 - \sqrt{2}\]

## PengRobinson78
[[PengRobinson78]]

- \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
- \[k = 0.37464 + 1.54226 \omega - 0.26992 \omega^2  \text{ where } \omega <=0.491\]
- \[k = 0.37464 + 1.48503 \omega - 0.16442 \omega^2  + 0.016666 \omega^3 \text{ where } \omega > 0.491\]
- \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
- \[b = 0.07779607r  R T_c/P_c\]
- \[\delta_1 = 1 + \sqrt{2}\]
- \[\delta_2 = 1 - \sqrt{2}\]

## RKPR
[[RKPR]]

The RKPR EoS extends the classical formulation of Cubic Equations of State by
freeing the parameter \(\delta_1\) and setting 
\(\delta_2 = \frac{1+\delta_1}{1-\delta_1}\).
This extra degree provides extra ways of implementing the equation in
comparison of other Cubic EoS (like PR and SRK) which are limited to
definition of their critical constants.

Besides that extra parameter, the RKRR includes another \(\alpha\)
function:

\[
 \alpha(T_r) = \left(\frac{3}{2+T_r}\right)^k
\]

In this implementation we take the simplest form which correlates
the extra parameter to the critical compressibility factor \(Z_c\) and
the \(k\) parameter of the \(\alpha\) function to \(Z_c\) and \(\omega\):

\[\delta_1 = d_1 + d_2 (d_3 - Z_c)^d_4 + d_5 (d_3 - Z_c) ^ d_6\]
\[k = (A_1  Z_c + A_0)\omega^2 + (B_1 Z_c + B_0)\omega + (C_1 Z_c + C_0)\]

It is also possible to include the parameters as optional arguments.