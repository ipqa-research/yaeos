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
\frac{d^2 \tau_{lk}}{dT^2} = \tau_{lk} \left(2 e_{lk} - \frac{c_{lk}}{T^{2}} +
\frac{2 b_{lk}}{T^{3}}\right)
$$

\(\qquad G^E\):

$$
\frac{\partial G^E}{\partial T} = \frac{G^E}{T} - RT \sum_k \frac{q_k
n_k \sum_l \theta_l \frac{\partial \tau_{lk}}{\partial T}}{\sum_l \theta_l
\tau_{lk}}
$$

$$
\frac{\partial G^E}{\partial T^2} = -R\left[T \sum_k q_k n_k
\left(\frac{(\sum_l \theta_l \frac{\partial^2 \tau_{lk}}{\partial T^2})}{\sum_l
\theta_l \tau_{lk}}
- \frac{(\sum_l \theta_l \frac{\partial \tau_{lk}}{\partial T})^2}{(\sum_l
- \theta_l \tau_{lk})^2}\right) + 2\left(\sum_k \frac{q_k n_k(\sum_l \theta_l
- \frac{\partial \tau_{lk}}{\partial T} )}{\sum_l \theta_l \tau_{lk}}\right)
- \right]
$$

## Cross temperature-compositional derivative

$$
\frac{\partial^2 G^E}{\partial n_i \partial T} = TODO
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
\begin{aligned}
& \frac{\partial \frac{G^E}{RT}}{\partial n_i} = \\

& + \frac{z}{2}{q}_{i}\ln{\left(\frac{\theta_{i}}{\phi_{i}} \right)} +
\frac{z}{2} \sum_k {n}_{k} {q}_{k} \frac{\left(\frac{d \theta_{k}}{d
{n}_{i}}\phi_{k} - \theta_{k} \frac{d \phi_{k}}{d {n}_{i}}\right)}{\theta_{k}
\phi_{k}} \\

& - {q}_{i} \ln{\left(\sum_l \theta_{l} {\tau}_{l,i} \right)} - \sum_k {n}_{k}
{q}_{k} \frac{\sum_l \frac{d \theta_{l}}{d {n}_{i}} {\tau}_{l,k}}{\sum_l
\theta_{l} {\tau}_{l,k}}

\end{aligned}
$$

\(\frac{\partial^2 \frac{G^E}{RT}}{\partial n_i \partial n_j} \) is obtained by
automatic differentiation.

## Examples


## References
1. Maurer, G., & Prausnitz, J. M. (1978). On the derivation and extension of
   the UNIQUAC equation. Fluid Phase Equilibria, 2(2), 91-99.
2. Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey. Chemical
   Thermodynamics for Process Simulation. 1st edition. Weinheim: Wiley-VCH,
   2012.
3. Caleb Bell and Contributors (2016-2024). Thermo: Chemical properties
   component of Chemical Engineering Design Library (ChEDL)
   https://github.com/CalebBell/thermo.