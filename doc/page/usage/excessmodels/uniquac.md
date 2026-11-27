---
title: UNIQUAC
---

[TOC]

# UNIQUAC

UNIQUAC (**uni**versal **qua**si**c**hemical) is an excess-Gibbs model for
nonideal liquid mixtures. It separates two reasons why real mixing differs from
ideal dilution: molecules occupy different amounts of space, and unlike
neighbors interact differently from like neighbors. Use it when you have
consistent molecular size/area parameters and fitted directional interaction
parameters for your components. It does not infer these from chemical names.

## Physical meaning and equation

The first two terms below are **combinatorial**: they account for molecular size
and shape. The final **residual** term describes energetic contacts. Here
"residual" is a part of the excess model, not an EOS departure from an ideal gas.
The complete \(G^E\) is relative to an ideal liquid solution; it is not the total
Gibbs energy of mixing. Differentiation gives \(\ln\gamma_i\), where activity
is \(x_i\gamma_i\). Fugacity also needs a reference fugacity, as explained in
the [shared introduction](index.html).

With the usual coordination number \(z=10\), UNIQUAC adds the three contributions:

$$ 
\frac{G^E}{RT} = \sum_k n_k \ln\frac{\phi_k}{x_k}
+ \frac{z}{2}\sum_k q_k n_k \ln\frac{\theta_k}{\phi_k}
- \sum_k q_k n_k \ln\left(\sum_l \theta_l \tau_{lk} \right)
$$

The mole fraction \(x_k\), volume fraction \(\phi_k\), and surface fraction
\(\theta_k\) use the same amounts but different weights:

$$x_k = \frac{n_k}{\sum_l n_l}$$

$$ \phi_k = \frac{r_k n_k}{\sum_l r_l n_l} $$

$$ \theta_k = \frac{q_k n_k}{\sum_l q_l n_l} $$

$$ \tau_{lk} = \exp \left[\frac{-\Delta U_{lk}}{R T} \right] $$

$$
\frac{-\Delta U_{lk}}{R T} = a_{lk}+\frac{b_{lk}}{T}+c_{lk}\ln T + d_{lk}T +
e_{lk}{T^2}
$$

The dimensionless molecular parameters \(r_k\) and \(q_k\) represent relative
volume and surface area, respectively. They are not measured volumes in L/mol
or areas in square metres. The contact factor \(\tau_{lk}\) weights how the
surface of component `l` interacts with component `k`; matrix direction matters.

Two conventions are worth stating explicitly, because both are easy to get
wrong when transcribing a published table.

- **Indices.** This is the *molecular* form of UNIQUAC: the sums run over
  components, not over functional groups. `r_k` and `q_k` belong to whole
  molecules, and `tau(l,k)` couples component `l` to component `k`. The
  interaction matrices are therefore `(nc,nc)`, with `nc` the number of
  components.
- **Reference state.** The residual term uses a pure-component reference, so
  the diagonal coefficients are normally zero: `tau(k,k) = 1` for a pure
  component, and the excess Gibbs energy of a pure component vanishes. Off-
  diagonal coefficients are directional and generally not symmetric.

## Constructor inputs and temperature convention

Use `setup_uniquac(qs, rs, aij, bij, cij, dij, eij)`, preferably with keywords.
The first array is **area** (`qs`) and the second is **volume** (`rs`). Both have
length `nc`, should be positive, and must use the same component order as `n`.
Each optional interaction array has shape `(nc,nc)`; omitted arrays are zero.

| Coefficient in `ln(tau)` | Units with `T` in kelvin |
|:--|:--|
| `aij` | Dimensionless |
| `bij/T` | `bij` in K |
| `cij*log(T)` | Dimensionless `cij`; logarithm uses the numerical kelvin value, equivalent to `ln(T/1 K)` |
| `dij*T` | `dij` in 1/K |
| `eij*T**2` | `eij` in 1/K² |

These are coefficients of **\(\ln\tau\)**, not of \(\Delta U\). In particular,
for a constant \(\Delta U_{ij}/R\) tabulated in kelvin, pass
`bij(i,j) = -Delta_U(i,j)/R`; the minus sign belongs in the conversion. Normally
all diagonal coefficients are zero so \(\tau_{ii}=1\). Off-diagonal coefficients
are directional and generally not symmetric; `tau(l,k)` enters the sum for
component `k`. Match your source's definition before deciding to transpose it.
Do not reuse NRTL coefficients simply because some arguments have similar names.

The runnable example is below the optional derivative reference. For a
structure-based model without molecular binary fits, see
[UNIFAC](unifac_models/index.html), which uses a related physical decomposition
but different parameters.

## Derivative reference

The methods evaluate derivatives for you. The rewritten expressions below avoid
ratios such as `0/0` in \(\phi_k/x_k\) when a component amount is zero, as long
as the total amount and weighted sums remain positive.

Some of the model's terms can be simplified to reduce the complexity of the
derivatives. This also allows the model to be evaluated at mole vectors `n`
where some of the amounts are zero.

$$\frac{\phi_k}{x_k} = \frac{n_T r_k}{\sum_l r_l n_l}$$

$$\frac{\theta_k}{\phi_k} = \frac{q_k \sum_l r_l n_l}{r_k \sum_l q_l n_l}$$

Being \(n_T \) the total number of moles in the system. The expression for the
Excess Gibbs free energy can be rewritten as:

$$ 
\frac{G^E}{RT} = \sum_k n_k \ln \left(\frac{n_T r_k}{\sum_l r_l n_l} \right)
+ \frac{z}{2}\sum_k q_k n_k \ln \left(\frac{q_k \sum_l r_l n_l}{r_k \sum_l q_l
n_l} \right) - \sum_k q_k n_k \ln\left(\sum_l \theta_l \tau_{lk} \right) 
$$


## Temperature derivatives

\(\qquad \tau_{lk}:\)

$$
\frac{d \tau_{lk}}{dT} = \tau_{lk} \left(2 T e_{lk} + d_{lk} + \frac{c_{lk}}{T}
- \frac{b_{lk}}{T^{2}}\right)
$$

$$
\frac{d^2 \tau_{lk}}{dT^2} = \tau_{lk} \left(2 T e_{lk} + d_{lk} + \frac{c_{lk}}{T}
- \frac{b_{lk}}{T^{2}}\right)^2 + \tau_{lk} \left(2 e_{lk} - \frac{c_{lk}}{T^{2}} +
\frac{2 b_{lk}}{T^{3}}\right)
$$

\(\qquad G^E\):

$$
\frac{\partial G^E}{\partial T} = \frac{G^E}{T} - RT \sum_k q_k n_k \frac{
\sum_l \theta_l \frac{\partial \tau_{lk}}{\partial T}}{\sum_l \theta_l
\tau_{lk}}
$$

$$
\frac{\partial^2 G^E}{\partial T^2} = -R\left[T \sum_k q_k n_k
\left(\frac{(\sum_l \theta_l \frac{\partial^2 \tau_{lk}}{\partial T^2})}{\sum_l
\theta_l \tau_{lk}}
- \frac{(\sum_l \theta_l \frac{\partial \tau_{lk}}{\partial T})^2}{(\sum_l
\theta_l \tau_{lk})^2}\right) + 2\left(\sum_k q_k n_k \frac{(\sum_l \theta_l
\frac{\partial \tau_{lk}}{\partial T} )}{\sum_l \theta_l \tau_{lk}}\right)
\right]
$$

## Cross temperature-compositional derivative

$$
\frac{\partial^2 G^E}{\partial n_i \partial T} = \frac{1}{T} \frac{\partial
G^E}{\partial n_i} - R T \left(q_i \frac{\sum_l \theta_l \frac{d \tau_{li}}{d
T}}{\sum_l \theta_l \tau_{li}} + \sum_k q_k n_k \left(\frac{(\sum_l \frac{d
\theta_l}{d n_i} \frac{d \tau_{lk}}{d T})(\sum_l \theta_l \tau_{lk}) - (\sum_l
\theta_l \frac{d \tau_{lk}}{d T})(\sum_l \frac{d \theta_l}{d n_i}
\tau_{lk})}{(\sum_l \theta_l \tau_{lk})^2} \right) \right)
$$


## Compositional derivatives

\(\qquad \phi_k\):

\(\phi_k \) derivatives are not really needed, but we also provide them.

$$
\frac{d \phi_k}{dn_i} = \begin{cases} - \frac{{n}_{i}
{r}_{i}^{2}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{2}} +
\frac{{r}_{i}}{\sum_l {n}_{l} {r}_{l}} & \text{for}\: i = k \\-
\frac{{n}_{k} {r}_{i} {r}_{k}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{2}}
& \text{otherwise} \end{cases}
$$

$$
\frac{d^2 \phi_k}{dn_i dn_j} = \begin{cases} \frac{2 {n}_{i}
{r}_{i}^{3}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{3}} - \frac{2
{r}_{i}^{2}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{2}} & \text{for}\: i
= k \wedge j = k \\\frac{2 {n}_{i} {r}_{i}^{2} {r}_{j}}{\left(\sum_l
{n}_{l} {r}_{l}\right)^{3}} - \frac{{r}_{i} {r}_{j}}{\left(\sum_l
{n}_{l} {r}_{l}\right)^{2}} & \text{for}\: i = k \wedge j \neq k \\\frac{2
{n}_{j} {r}_{i} {r}_{j}^{2}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{3}} -
\frac{{r}_{i} {r}_{j}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{2}} &
\text{for}\: j = k \wedge i \neq k \\\frac{2 {n}_{k} {r}_{i} {r}_{j}
{r}_{k}}{\left(\sum_l {n}_{l} {r}_{l}\right)^{3}} & \text{otherwise}
\end{cases}
$$

\(\qquad \theta_k\):

$$
\frac{d \theta_k}{dn_i} = \begin{cases} - \frac{{n}_{i}
{q}_{i}^{2}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{2}} +
\frac{{q}_{i}}{\sum_l {n}_{l} {q}_{l}} & \text{for}\: i = k \\-
\frac{{n}_{k} {q}_{i} {q}_{k}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{2}}
& \text{otherwise} \end{cases}
$$

$$
\frac{d^2 \theta_k}{dn_i dn_j} = \begin{cases} \frac{2 {n}_{i}
{q}_{i}^{3}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{3}} - \frac{2
{q}_{i}^{2}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{2}} & \text{for}\: i
= k \wedge j = k \\\frac{2 {n}_{i} {q}_{i}^{2} {q}_{j}}{\left(\sum_l
{n}_{l} {q}_{l}\right)^{3}} - \frac{{q}_{i} {q}_{j}}{\left(\sum_l
{n}_{l} {q}_{l}\right)^{2}} & \text{for}\: i = k \wedge j \neq k \\\frac{2
{n}_{j} {q}_{i} {q}_{j}^{2}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{3}} -
\frac{{q}_{i} {q}_{j}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{2}} &
\text{for}\: j = k \wedge i \neq k \\\frac{2 {n}_{k} {q}_{i} {q}_{j}
{q}_{k}}{\left(\sum_l {n}_{l} {q}_{l}\right)^{3}} & \text{otherwise}
\end{cases}
$$

\(\qquad G^E\):

$$
\frac{\partial \frac{G^E}{RT}}{\partial n_i} = \ln 
\left(\frac{n_T r_i}{\sum_l r_l n_l} \right) + \sum_k n_k 
\left(\frac{\sum_l r_l n_l - n_T r_i}{n_T \sum_l r_l n_l}\right) +
\frac{z}{2} q_i \ln \left(\frac{q_i \sum_l r_l n_l}{r_i \sum_l q_l n_l} 
\right) + \frac{z}{2} \sum_k n_k q_k \left(\frac{r_i \sum_l n_l q_l - q_i
\sum_l n_l r_l}{(\sum_l n_l q_l) (\sum_l n_l r_l)} \right) - {q}_{i}
\ln{\left(\sum_l \theta_{l} {\tau}_{li} \right)} - \sum_k {n}_{k} {q}_{k}
\frac{\sum_l \frac{d \theta_{l}}{d {n}_{i}} {\tau}_{lk}}{\sum_l \theta_{l}
{\tau}_{lk}}
$$

Differentiating each term of the first compositional derivative respect to
\(n_j\) we get:

\(\frac{\partial^2 \frac{G^E}{RT}}{\partial n_i \partial n_j} =\)

$$
\frac{\sum_l n_l r_l - n_T r_j}{n_T \sum_l n_l r_l} 
$$

$$
+ \frac{\sum_l n_l r_l - n_T r_i}{n_T \sum_l r_l n_l}
+ \sum_k n_k \left(\frac{r_i r_j}{(\sum_l n_l r_l)^2} - \frac{1}{n_T^2} \right)
$$

$$
+ \frac{z}{2} q_i \left( \frac{r_j \sum_l n_l q_l - q_j \sum_l n_l r_l}
{(\sum_l n_l q_l)(\sum_l n_l r_l)}\right)
$$

$$
+ \frac{z}{2} q_j \left( \frac{r_i \sum_l n_l q_l - q_i \sum_l n_l r_l}
{(\sum_l n_l q_l)(\sum_l n_l r_l)}\right) + \frac{z}{2} \sum_k n_k q_k \left(
\frac{q_i q_j}{(\sum_l n_l q_l)^2} - \frac{r_i r_j}{(\sum_l n_l r_l)^2} \right)
$$

$$
- q_i \left( \frac{\sum_l \frac{d \theta_l}{d n_j} \tau_{li}}{\sum_l \theta_l 
\tau_{li}} \right)
$$

$$
- {q}_{j} \frac{\sum_l \frac{d \theta_{l}}{d {n}_{i}} {\tau}_{lj}}{\sum_l
\theta_{l}{\tau}_{lj}} - \sum_k {n}_{k} {q}_{k} \frac{\left(\sum_l
\frac{d^2\theta_l}{dn_idn_j} \tau_{lk} \right) \left(\sum_l
\theta_l \tau_{lk} \right) - \left(\sum_l \frac{d\theta_l}{dn_i}
\tau_{lk} \right) \left(\sum_l \frac{d\theta_l}{dn_j} \tau_{lk}
\right)}{(\sum_l \theta_{l} {\tau}_{lk})^2}
$$

## Examples
Example from: Gmehling et al. (2012) [2]

An example of a mixture of water–ethanol–benzene at 298.15 K with
constant \(\frac{\Delta U}{R}\) [K]:

|Water|Ethanol|Benzene|
|---------|--------|---------|
| 0       | 526.02 | 309.64  |
| −318.06 | 0      | −91.532 |
| 1325.1  | 302.57 | 0       |

Read that table as \(\Delta U_{lk}/R\). Because the model expects the
coefficients of \(\ln\tau\), and \(\ln\tau_{lk}=-\Delta U_{lk}/(RT)\), a constant
energy difference becomes a coefficient of \(1/T\), namely
`bij(l,k) = -Delta_U(l,k)/R`. That is exactly the minus sign in the code
below. It is the single most common transcription error: entering the table
values directly produces a model with the wrong sign of interaction.

```fortran
program uniquac_activity
use yaeos, only: pr, setup_uniquac, UNIQUAC
implicit none

integer, parameter :: nc = 3

real(pr) :: rs(nc)
real(pr) :: qs(nc)
real(pr) :: b(nc, nc)
real(pr) :: n(nc)

real(pr) :: ln_gammas(nc)
real(pr) :: T
real(pr) :: Ge

type(UNIQUAC) :: model

rs = [0.92_pr, 2.1055_pr, 3.1878_pr]
qs = [1.4_pr, 1.972_pr, 2.4_pr]

T = 298.15_pr

! Calculate bij from DUij. We need -DU/R to get bij
b(1,:) = [0.0_pr, -526.02_pr, -309.64_pr]
b(2,:) = [318.06_pr, 0.0_pr, 91.532_pr]
b(3,:) = [-1325.1_pr, -302.57_pr, 0.0_pr]

model = setup_uniquac(qs=qs, rs=rs, bij=b)

n = [2.0_pr, 2.0_pr, 8.0_pr]

call model%ln_activity_coefficient(n, T, ln_gammas)

call model%excess_gibbs(n, T, Ge=Ge)
print *, "gamma: water, ethanol, benzene = ", exp(ln_gammas)
print *, "GE per mole [bar L/mol] = ", Ge/sum(n)
print *, "Composition = ", n/sum(n)
end program uniquac_activity

```

This is a complete program. Table rows and explicit `b(i,:)` assignments use the
same water–ethanol–benzene order. Explicit rows avoid the column-major ordering
trap of `reshape`. The amounts total 12 mol, so the composition is not the raw
`n` vector. Coefficients are dimensionless; they may be above or below 1, and
the output of `ln_activity_coefficient` must be exponentiated before interpreting
it as `gamma`. No claim of a stable, single liquid follows from this property
evaluation of a trial composition.

The printed coefficients show how differently the components respond: water's
coefficient is far above 1 while ethanol's is below 1, and benzene sits near
the middle. Comparing these numbers with the values reported for the same
system in the source reference is the most direct check that your parameter
table and sign convention were transcribed correctly. A sign error in `bij`
reverses this pattern, so it shows up immediately here.

### Adding temperature-dependent parameters

The five optional matrices allow a richer temperature dependence. Passing
`dij` or `eij` adds a linear or quadratic term to \(\ln\tau\), and `cij` adds a
logarithmic term. For example, to give one pair a temperature-dependent
interaction, add to the program above:

```fortran
real(pr) :: dij(nc, nc)
real(pr) :: eij(nc, nc)

dij = 0.0_pr
eij = 0.0_pr
dij(1, 2) = 1.0e-3_pr

model = setup_uniquac(qs=qs, rs=rs, bij=b, dij=dij, eij=eij)
```

The units follow the table in [Constructor inputs](#constructor-inputs-and-temperature-convention):
`dij` in 1/K and `eij` in 1/K². `aij` adds a constant offset to \(\ln\tau\).
Use these terms only when your parameter source provides them; they are not
defaults that improve accuracy by themselves.

## Checks, common errors, and limits

The most consequential mistake with UNIQUAC is a unit or sign error in the
interaction parameters, because the model still returns a plausible-looking
number. Use kelvin, never Celsius. If your energy table is in J/mol, convert
it with a matching gas constant before setting `bij` in kelvin, and remember the
minus sign described earlier. Do not swap `qs` and `rs`: changing molecular
geometry is not equivalent to changing binary interaction energies, and the
two errors are easy to make because the argument order is area first.

A useful sanity check is the ideal limit. Unlike NRTL, setting the
interaction matrices to zero does **not** generally make UNIQUAC ideal, because
the combinatorial contribution survives: different-sized molecules still
occupy different amounts of space. Identical `r` and `q` for every component
together with zero interaction coefficients does give the ideal result, and is
worth checking once when you first set up a model.

Two further checks apply to every variant. Scaling all amounts must leave the
activity coefficients unchanged. At a pure composition, the present
component's coefficient approaches 1 for properly referenced parameters;
coefficients of absent solutes need not approach 1, because a solute at
infinite dilution has its own limiting value.

Finally, a fitted temperature dependence is not a license for extrapolation.
Very large exponents can overflow, and the model has no general electrolyte or
near-critical treatment and supplies no density or vapor properties.

Next, evaluate [excess enthalpy and other shared properties](index.html) and
compare with your parameter source. For phase splits, use
[stability and equilibrium methods](../phase_equilibrium/index.html), not the
sign of `Ge` or the size of a single coefficient alone.


## References
1. Maurer, G., & Prausnitz, J. M. (1978). On the derivation and extension of
   the UNIQUAC equation. Fluid Phase Equilibria, 2(2), 91-99.
2. Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey. Chemical
   Thermodynamics for Process Simulation. 1st edition. Weinheim: Wiley-VCH,
   2012.
3. Caleb Bell and Contributors (2016-2024). Thermo: Chemical properties
   component of Chemical Engineering Design Library (ChEDL)
   https://github.com/CalebBell/thermo.
