---
title: NRTL-HV
---

# NRTL-HV


The NRLT model presented by Huron and Vidal is a variation of the original
non-random two-liquid model. This variation is intended to be used combined with
cubic Equations of State, using the Huron-Vidal mixing rule.  The major
difference with respect to the original NRTL model is that it includes the
repulsive parameter of the cubic EoS in the expression for the excess Gibbs
energy. This allows the model to be simplified to the classic Van der Waals
1-fluid quadratic mixing rules when giving special values to the parameters.

The general expression for the excess Gibbs energy is given by:

\begin{equation}
    \frac{G^E}{RT} = \sum\limits_i n_i 
    \frac
    {\sum\limits_j n_j b_j E_{j,i} \tau_{j,i}}
    {\sum\limits_j n_j b_j E_{j, i}}
\end{equation}

Where the parameters are defined as follows:

\begin{align}
    \tau_{j, i}(T) & = \frac{g_{j,i} - g_{i,i}}{RT} \\
    E_{j, i} &= \exp\left(-\alpha_{j, i} \tau_{j, i} \right)
\end{align}

In this implementation:

\begin{equation}
g_{j,i} = g_{j,i}^0 + g_{j,i}^T \cdot T
\end{equation}

Where \(g_{j,i}\) is the interaction parameter between components \(j\) and \(i\),
which can be a function of temperature, and \(\alpha_{j,i}\) is the non-randomness
parameters which is usually symmetric with values between 0 and 0.5.

## Derivatives

In this section we present the first and second order derivatives of the excess
Gibbs energy. To simplify the expressions we will define the following
parameters:

\begin{align}
    \xi_{ji} & = E_{ji} \tau_{ji} b_j n_j \\
    \theta{ji} & = {E_{ji} b_j n_j} \\
    \Omega{ji} &= E_{ji} b_j \\
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
    \frac{1}{RT}\frac{dG^E}{dn_i} = \left(\sum_{i=1}^{NC} \frac{{E}_{k,i}
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
    \frac{1}{RT}\frac{dG^E}{dn_{ij}} = 
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
    \frac{\partial G^E}{\partial n_i \partial T}  = & \\
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


