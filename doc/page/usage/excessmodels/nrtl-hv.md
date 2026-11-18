---
title: NRTL-HV
---

[TOC]

# NRTL-HV


The NRTL model modified by Huron and Vidal is a variation of the original
non-random two-liquid model. This variation is intended to be used combined with
cubic Equations of State, using the Huron-Vidal mixing rule.  The major
difference with respect to the original NRTL model is that it includes the
repulsive parameter of the cubic EoS in the expression for the excess Gibbs
energy. This allows the model to be simplified to the classic Van der Waals
1-fluid quadratic mixing rules when giving special values to the parameters.

## When to use it

Choose this model when your parameter set was developed for a Huron-Vidal
cubic-EOS formulation. Ordinary [NRTL](nrtl.html) weights neighbors by their
amounts and interaction factors; NRTL-HV additionally weights them by pure-fluid
**covolumes** \(b_j\), which represent excluded-volume size in the cubic EOS.
This is not a liquid molar-volume calculation. In `yaeos` the public type is
[[NRTLHV]], while the ordinary type is [[NRTL]].

NRTL-HV itself remains a [[GeModel]]: its derivatives give activity coefficients
relative to an ideal liquid, not fugacity coefficients relative to an ideal gas.
An activity is \(x_i\gamma_i\); a liquid fugacity also requires its pure-liquid
reference fugacity. Coupling the excess model to a cubic EOS is a separate step
through the [mixing rule](../eos/cubics/mixing.html), not a side effect of calling
the constructor.

## Equation and actual parameter convention

The general expression for the excess Gibbs energy is given by:

\begin{equation}
    \frac{G^E}{RT} = \sum\limits_i n_i 
    \frac
    {\sum\limits_j n_j b_j E_{j,i} \tau_{j,i}}
    {\sum\limits_j n_j b_j E_{j, i}}
\end{equation}

Where the parameters are defined as follows:

\begin{align}
    \tau_{j, i}(T) & = \frac{\Delta g_{j,i}(T)}{RT} \\
    E_{j, i} &= \exp\left(-\alpha_{j, i} \tau_{j, i} \right)
\end{align}

In this implementation:

\begin{equation}
\Delta g_{j,i}(T) = g_{j,i}^0 + g_{j,i}^T \cdot T
\end{equation}

Here \(\Delta g_{j,i}\) is the **already-referenced** interaction difference.
If a source supplies raw energies \(g_{j,i}-g_{i,i}\), subtract the reference
before constructing this model: the implementation does **not** subtract a
diagonal energy. The non-randomness parameter \(\alpha_{j,i}\) controls the local
weighting and is often symmetric, but must match the fitted parameter set.

| Constructor argument | Meaning and units |
|:--|:--|
| `b(:)` | Pure-component cubic covolumes, L/mol |
| `alpha(:,:)` | Dimensionless non-randomness parameters |
| `gji0(:,:)` | Constant referenced energy differences, bar L/mol |
| `gjiT(:,:)` | Coefficients multiplying `T`, bar L/(mol K) |

Thus \(\tau_{ji}=\texttt{gjiT(j,i)}/R+\texttt{gji0(j,i)}/(RT)\), with `T`
in kelvin and the library gas constant `R` in bar L/(mol K). Divide J/mol by 100
to obtain bar L/mol (and likewise for J/(mol K)). These arrays are **not** the
dimensionless `a` and kelvin-valued `b` matrices of the ordinary NRTL constructor.
All matrices have shape `(nc,nc)` and need not be symmetric. `(:,i)` describes
contacts in the local average around component `i`. Use zero diagonal energy
differences and keep component order consistent with `n` and the covolume vector.

## Complete property example

The following artificial binary example tests the API and its ideal-interaction
limit; it is not a parameter set for named chemicals. For EOS work, obtain `b`
from the same pure-component cubic model used for coupling.

```fortran
program nrtlhv_activity
use yaeos, only: pr, NRTLHV
implicit none

type(NRTLHV) :: model
real(pr) :: b(2)
real(pr) :: alpha(2, 2)
real(pr) :: gji0(2, 2)
real(pr) :: gjiT(2, 2)
real(pr) :: n(2)
real(pr) :: T
real(pr) :: ln_gamma(2)
real(pr) :: Ge

b = [0.04_pr, 0.06_pr]
alpha = 0.3_pr
gji0 = 0.0_pr
gji0(1, 2) = 10.0_pr
gji0(2, 1) = -5.0_pr
gjiT = 0.0_pr
n = [0.4_pr, 0.6_pr]
T = 300.0_pr

model = NRTLHV(b=b, alpha=alpha, gji0=gji0, gjiT=gjiT)
call model%ln_activity_coefficient(n, T, ln_gamma)
call model%excess_gibbs(n, T, Ge=Ge)
print *, "gamma, components 1 and 2 = ", exp(ln_gamma)
print *, "GE [bar L] = ", Ge

gji0 = 0.0_pr
model = NRTLHV(b=b, alpha=alpha, gji0=gji0, gjiT=gjiT)
call model%ln_activity_coefficient(n, T, ln_gamma)
print *, "Zero-interaction ln(gamma), near zero = ", ln_gamma
end program nrtlhv_activity
```

The output coefficient order follows `n`. Coefficients must be positive, not
necessarily less than 1. With both energy matrices zero, \(G^E=0\) and the
coefficients are 1, even when the covolumes differ. This ideal-*excess* check is
not the same as the special parameter choices needed to recover quadratic
mixing in an EOS.

### Cross-check against ordinary NRTL

When every component has the same covolume, the covolume weighting drops out
of the local averages and NRTL-HV reduces to ordinary NRTL for the same
\(\tau\) and \(\alpha\). Because the two models use different parameter forms,
the conversion is:

| Ordinary NRTL | NRTL-HV equivalent |
|:--|:--|
| `a(i,j)` (dimensionless) | `gjiT(i,j) = a(i,j) * R` |
| `b(i,j)` (K) | `gji0(i,j) = b(i,j) * R` |
| `c(i,j)` (dimensionless) | `alpha(i,j) = c(i,j)` |
| — | `b(i) = 1` for every component |

With those substitutions the two models give identical activity coefficients.
This is a useful diagnostic when converting a parameter set: if the results
differ, a unit or a matrix direction is wrong.

### Checks, limits, and next steps

Use positive covolumes and nonnegative mole amounts with a positive total.
Wrong energy units, reversed matrix directions, and subtracting the diagonal
twice are particularly damaging mistakes. Equal covolumes reduce the weighting
to ordinary NRTL for the same \(\tau\) and \(\alpha\), which offers another
cross-check. Do not transplant a fitted ordinary-NRTL parameter set into an
unequal-covolume model and assume unchanged predictions.

The special parameter choices that recover the classic van der Waals one-fluid
quadratic mixing rule inside a cubic EOS (a specific set of diagonal and
off-diagonal energies, and zero non-randomness) are implemented in the
[Huron-Vidal mixing rule](../eos/cubics/mixing.html), not by this model object.

No density, vapor phase, ionic chemistry, or pressure-dependent excess properties
are supplied by this standalone model. Validate the fitted temperature range,
then continue to [shared property methods](index.html) and
[cubic mixing rules](../eos/cubics/mixing.html). The derivative details below are
optional background; application code should call the model methods.

## Derivatives

In this section we present the first and second order derivatives of the excess
Gibbs energy. To simplify the expressions we will define the following
parameters:

\begin{align}
    \xi_{ji} & = E_{ji} \tau_{ji} b_j n_j \\
    \theta_{ji} & = {E_{ji} b_j n_j} \\
    \Omega_{ji} &= E_{ji} b_j \\
    \eta_{ji} &= E_{ji} \tau_{ji} b_j
\end{align}

With this definitions, the excess Gibbs energy can be rewritten as:

\begin{equation}
    \frac{G^E}{RT} = \sum\limits_i n_i \frac
    {\sum\limits_j \xi_{j,i}}{\sum\limits_j \theta_{j,i}}
\end{equation}

### Compositional derivatives

The first order derivatives of the excess Gibbs energy with respect to the
number of moles of component \(i\) are given by:

\begin{equation}
    \frac{1}{RT}\frac{\partial G^E}{\partial n_k} = \left(\sum_{i=1}^{NC} \frac{{E}_{k,i}
    {\tau}_{k,i} {b}_{k} {n}_{i}}{\sum_{j=1}^{NC} {E}_{j,i} {b}_{j} {n}_{j}} +
    \sum_{i=1}^{NC} - \frac{{E}_{k,i} {b}_{k} {n}_{i} \sum_{j=1}^{NC} {E}_{j,i}
    {\tau}_{j,i} {b}_{j} {n}_{j}}{\left(\sum_{j=1}^{NC} {E}_{j,i} {b}_{j}
    {n}_{j}\right)^{2}} + \frac{\sum_{j=1}^{NC} {E}_{j,k} {\tau}_{j,k} {b}_{j}
    {n}_{j}}{\sum_{j=1}^{NC} {E}_{j,k} {b}_{j} {n}_{j}}\right)
\end{equation}

Which can be rewritten in terms of the defined parameters as:

\begin{equation}
\frac{1}{RT}\frac{dG^E}{dn_i} =
\frac{\sum\limits_l \xi_{l, i}}{\sum\limits_l \theta_{l, i}} 
+ \sum\limits_k n_k 
\left(
    \frac{\eta_{i, k}}{\sum\limits_l \theta_{l, k}} - \frac{\Omega_{i, k}
    \sum\limits_l \xi_{l, k}}{\left(\sum\limits_l \theta_{l,k}\right)^2}
\right)
\end{equation}

The second order derivatives of the excess Gibbs energy with respect to the number of moles of components \(i\) and \(j\) are given by:

\begin{equation}
    \frac{1}{RT}\frac{\partial^2 G^E}{\partial n_i\partial n_j} =
    \begin{aligned}
    - \frac{\Omega_{ji} \sum\limits_k{\xi_{ki}}}{\left(\sum\limits_k \theta_{ki}\right)^2}
    - \frac{\Omega_{ij} \sum\limits_k \xi_{kj}}{\left(\sum\limits_k \theta_{kj}\right)^2}
    + \frac{\eta_{ji}}{\sum\limits_k \theta_{ki}}
    + \frac{\eta_{ij}}{\sum\limits_k \theta_{kj}}
    + \sum\limits_k 
    \left(
        \frac{2 n_k \Omega_{ik} \Omega_{jk} \sum\limits_l \xi_{lk}}{\left(\sum\limits_l \theta_{lk}\right)^3}
        - \frac{n_k \Omega_{ik} \eta_{jk}}{\left(\sum\limits_l \theta_{lk}\right)^2}
        - \frac{n_k \Omega_{jk} \eta_{ik}}{\left(\sum\limits_l \theta_{lk}\right)^2}
    \right)
    \end{aligned}
\end{equation}

### Temperature derivatives

To simplify the equations, we present the derivatives of each variable as
\(x^T\) for first order derivatives wrt to temperature and \(x^{TT}\) for
second order derivatives.

First derivative
\begin{equation}
    \frac{\partial G^E}{\partial T} = \frac{G^E}{T} + RT \sum\limits_i n_i \left( \frac{\sum\limits_j \xi^T_{j,i}}{\sum\limits_j \theta_{j, i}} - \frac{\sum\limits_j \xi_{j,i} \sum\limits_j \theta^T_{j,i}}{\left(\sum\limits_j{\theta_{j,i}}\right)^2}\right)
\end{equation}

Second derivative

\begin{equation}
    \begin{aligned}
    \frac{\partial^2 G^E}{\partial T^2} =& 
    \left( -\frac{G^E}{T^2} + \frac{\frac{\partial G^E}{\partial T}}{T} \right) 
    + R \sum\limits_i n_i \left( \frac{\sum\limits_j \xi^T_{j,i}}{\sum\limits_j \theta_{j, i}} - \frac{\sum\limits_j \xi_{j,i} \sum\limits_j \theta^T_{j,i}}{\left(\sum\limits_j{\theta_{j,i}}\right)^2}\right) \\
    &+ RT \sum\limits_i n_i 
    \left(
        \frac{\sum\limits_j \xi_{j,i}^{TT}}{\sum\limits_j \theta_{j,i}} - \frac{\sum\limits_j \xi^{T}_{j,i} \sum\limits_j \theta^T_{j,i}}{\left(\sum\limits_j{\theta_{j,i}}\right)^2}
        - \left(
            \frac{\sum\limits_j \xi^T_{j,i}\sum\limits_j\theta^T_{j,i} + \sum\limits_j\xi_{j,i}\sum\limits_j\theta^{TT}_{j,i}}{\left(\sum\limits_j \theta_{j,i}\right)^2} \right. \right. \\
        & \left. \left. 
            - 2 \frac{\left(\sum\limits_j\theta^T_{j,i}\right)^2 \sum\limits_j\xi_{j,i}}{\left(\sum\limits_j \theta_{j,i}\right)^3}
        \right)
        \right)
    \end{aligned}
\end{equation}

### Crossed derivative

\begin{equation}
\begin{aligned}
    \frac{\partial^2 G^E}{\partial n_i \partial T}  = & \\
    & \frac{\frac{\partial G^E}{\partial n_i}}{T}  + RT \left( 
       \frac{\sum\limits_l \xi_{l,i}^T}{\sum\limits_l\theta_{l,i}} - \sum\limits_l \xi_{l, i} \frac{\sum\limits_l \theta^T_{l,i}}{\left(\sum\limits_l\theta_{l,i}\right)^2} \right.\\
    & + \sum\limits_k^{NC} n_k 
    \left. \left[
        \frac{\eta_{i,k}^T}{\sum\limits_l \theta_{l,k}} - \eta_{i,k}\frac{\sum\limits_l \theta^T_{l,k}}{\left(\sum\limits_l \theta_{l,k}\right)^2} 
        - \left(
            \frac{\Omega_{i,k} \sum\limits_l \xi_{l,k}^T + \Omega_{i, k}^T \sum\limits_l\xi_{l,k}}{\left(\sum\limits_l \theta_{l,k}\right)^2} 
            - 2 \frac{\sum\limits_l \theta_{l,k}^T \Omega_{i,k} \sum\limits_l \xi_{l,k}}{\left(\sum\limits_l \theta_{l,k}\right)^3}
        \right)
    \right]
    \right)
\end{aligned}
\end{equation}

## References

1. Huron, M.-J., and Vidal, J. (1979). New mixing rules in simple equations of
   state for representing vapour-liquid equilibria of strongly non-ideal
   mixtures. Fluid Phase Equilibria, 3, 255–271.
