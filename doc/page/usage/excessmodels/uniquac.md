---
title: UNIQUAC
---

[TOC]

# UNIQUAC

UNIQUAC (**uni**versal **qua**si**c**hemical) Excess Gibbs free energy model.

$$ 
\frac{G^E}{RT} = \sum_k n_k \ln\frac{\phi_k}{x_k}
+ \frac{z}{2}\sum_k q_k n_k \ln\frac{\theta_k}{\phi_k}
- \sum_k q_k n_k \ln\left(\sum_l \theta_l \tau_{lk} \right)
$$

With:

$$x_k = \frac{n_k}{\sum_l n_l}$$

$$ \phi_k = \frac{r_k n_k}{\sum_l r_l n_l} $$

$$ \theta_k = \frac{q_k n_k}{\sum_l q_l n_l} $$

$$ \tau_{lk} = \exp \left[\frac{-\Delta U_{lk}}{R T} \right] $$

$$
\frac{-\Delta U_{lk}}{R T} = a_{lk}+\frac{b_{lk}}{T}+c_{lk}\ln T + d_{lk}T +
e_{lk}{T^2}
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
\frac{\partial G^E}{\partial T^2} = -R\left[T \sum_k q_k n_k
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
\frac{\partial \frac{G^E}{RT}}{\partial n_i} = \ln \left(\frac{\phi_i}{x_i}
\right) + \sum_k n_k \left(\frac{\frac{d \phi_k}{dn_i}}{\phi_k} -
\frac{\frac{dx_k}{dn_i}}{x_k}\right) +
\frac{z}{2}{q}_{i}\ln{\left(\frac{\theta_{i}}{\phi_{i}} \right)} + \frac{z}{2}
\sum_k {n}_{k} {q}_{k} \left(\frac{\frac{d \theta_{k}}{d {n}_{i}}}{\theta_k} -
\frac{\frac{d \phi_{k}}{d {n}_{i}}}{\phi_k} \right) - {q}_{i}
\ln{\left(\sum_l \theta_{l} {\tau}_{li} \right)} - \sum_k {n}_{k} {q}_{k}
\frac{\sum_l \frac{d \theta_{l}}{d {n}_{i}} {\tau}_{lk}}{\sum_l \theta_{l}
{\tau}_{lk}}
$$


Differentiating each term of the first compositional derivative respect to
\(n_j\) we get:

\(\frac{\partial \frac{G^E}{RT}}{\partial n_i \partial n_j} =\)

$$
\frac{\frac{d \phi_i}{d n_j}}{\phi_i} - \frac{\frac{d x_i}{d n_j}}{x_i} 
$$

$$
+\frac{\frac{d \phi_j}{dn_i}}{\phi_j} -
\frac{\frac{dx_j}{dn_i}}{x_j}
+ \sum_k n_k \left(\frac{\frac{d^2\phi_k}{dn_i dn_j} \phi_k -
  \frac{d\phi_k}{dn_i} \frac{d\phi_k}{dn_j}}{\phi_k^2} \right)
- \sum_k n_k \left(\frac{\frac{d^2x_k}{dn_i dn_j} x_k -
  \frac{dx_k}{dn_i} \frac{dx_k}{dn_j}}{x_k^2} \right)
$$

$$
+ \frac{z}{2} q_i \left( \frac{\frac{d \theta_i}{d n_j}}{\theta_i} - \frac
{\frac{d \phi_i}{d n_j}}{\phi_i} \right)
$$

$$
+ \frac{z}{2} q_j \left( \frac{\frac{d \theta_j}{d n_i}}{\theta_j} - \frac
{\frac{d \phi_j}{d n_i}}{\phi_j} \right)
+ \frac{z}{2} \sum_k n_k q_k \left(\frac{\frac{d^2\theta_k}{dn_i dn_j} \theta_k -
  \frac{d\theta_k}{dn_i} \frac{d\theta_k}{dn_j}}{\theta_k^2} \right)
- \frac{z}{2} \sum_k n_k q_k \left(\frac{\frac{d^2\phi_k}{dn_i dn_j} \phi_k -
  \frac{d\phi_k}{dn_i} \frac{d\phi_k}{dn_j}}{\phi_k^2} \right)
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

An example of having a mixture of Water-Ethanol-Bezene at 298.15 K with 
constant \(\frac{\Delta U}{R}\) [K]:

|Water|Ethanol|Benzene|
|---------|--------|---------|
| 0       | 526.02 | 309.64  |
| −318.06 | 0      | −91.532 |
| 1325.1  | 302.57 | 0       |


```fortran
use yaeos, only: pr, setup_uniquac, UNIQUAC

integer, parameter :: nc = 3

real(pr) :: rs(nc), qs(nc)
real(pr) :: b(nc, nc)
real(pr) :: n(nc)

real(pr) :: ln_gammas(nc), T

type(UNIQUAC) :: model

rs = [0.92_pr, 2.1055_pr, 3.1878_pr]
qs = [1.4_pr, 1.972_pr, 2.4_pr]

T = 298.15_pr

! Calculate bij from DUij. We need -DU/R to get bij
b(1,:) = [0.0_pr, -526.02_pr, -309.64_pr]
b(2,:) = [318.06_pr, 0.0_pr, 91.532_pr]
b(3,:) = [-1325.1_pr, -302.57_pr, 0.0_pr]

model = setup_uniquac(qs, rs, bij=b)

n = [2.0_pr, 2.0_pr, 8.0_pr]

call model%ln_activity_coefficient(n, T, ln_gammas)

print *, exp(ln_gammas) ! [8.856, 0.860, 1.425]

```


## References
1. Maurer, G., & Prausnitz, J. M. (1978). On the derivation and extension of
   the UNIQUAC equation. Fluid Phase Equilibria, 2(2), 91-99.
2. Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey. Chemical
   Thermodynamics for Process Simulation. 1st edition. Weinheim: Wiley-VCH,
   2012.
3. Caleb Bell and Contributors (2016-2024). Thermo: Chemical properties
   component of Chemical Engineering Design Library (ChEDL)
   https://github.com/CalebBell/thermo.