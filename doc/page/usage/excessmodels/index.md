---
title: Gibbs Excess Models
---

[TOC]

# Gibbs Excess Models

Excess properties are properties of mixtures which quantify the non-ideal
behavior of real mixtures. They are defined as the difference between the value
of the property in a real mixture and the value that would exist in an ideal
solution under the same conditions.

$$
    G^E = G - G^{I,S} = \Delta G_{mix}
$$

$$
    H^E = H - H^{I,S} = \Delta H_{mix}
$$

$$
    S^E = H - H^{I,S} = \Delta S_{mix}
$$

Gibbs excess models are defined in terms of a mathematical expression that
describes the excess Gibbs energy of a mixture as a function of the composition
of the components in the mixture and temperature. Then, the other excess
properties can be derived from the excess Gibbs energy derivatives.


## Routines and Methods
### \( \ln \gamma_i \)

Activity coefficients are calculated from the excess Gibbs energy using the
following expression:

$$
    \ln \gamma_i = \frac{\partial \left(\frac{G^E}{RT} \right)}{\partial n_i}
$$

#### \( \ln \gamma_i \) derivatives

The derivatives of the activity coefficients with respect to the mole numbers:

$$
    \frac{\partial \ln \gamma_i}{\partial n_j} = \frac{\partial^2 \left(\frac{G^E}{RT} \right)}{\partial n_i \partial n_j}
$$

The temperature derivative of the activity coefficients:

$$
    \frac{\partial \ln \gamma_i}{\partial T} = \frac{\partial^2 
    \left(\frac{G^E}{RT} \right)}{\partial n_i \partial T} = 
    \frac{1}{RT} \left(
    \frac{\partial^2 G^E}{\partial n_i \partial T} - 
    \frac{1}{T} \frac{\partial G^E}{\partial n_i} \right) 
$$

You may notice that `yaeos` calculates the derivatives of the \(\ln \gamma \).
If for some reason you need the derivatives of the activity coefficients, you
may find useful the following expressions:

$$
    \gamma_i = e^{\ln \gamma_i}
$$

$$
    \frac{\partial \gamma_i}{\partial n_j} = \frac{\partial \ln \gamma_i}{\partial n_j} e^{\ln \gamma_i}
$$

$$
    \frac{\partial \gamma_i}{\partial T} = 
    \frac{\partial \ln \gamma_i}{\partial T} e^{\ln \gamma_i}
$$


#### Example

```fortran
program main
use yaeos, only: pr
use yaeos, only: Groups, setup_unifac, UNIFAC

type(UNIFAC) :: model

integer, parameter :: nc = 3

type(Groups) :: molecules(nc)

real(pr) :: n(nc), T, ln_gamma(nc), dln_gamma_dT(nc)
real(pr) :: dln_gamma_dn(nc, nc)

T = 298.15_pr
n = [20.0_pr, 70.0_pr, 10.0_pr]

! ! Ethane [CH3]
molecules(1)%groups_ids = [1]
molecules(1)%number_of_groups = [2]

! ! Ethanol [CH3, CH2, OH]
molecules(2)%groups_ids = [1, 2, 14]
molecules(2)%number_of_groups = [1, 1, 1]

! ! Methylamine [CH3-NH2]
molecules(3)%groups_ids = [28]
molecules(3)%number_of_groups = [1]

! setup UNIFAC model
model = setup_unifac(molecules)

! Calculate ln_gamma and its derivatives
call model%ln_gamma(n, T, ln_gamma, dln_gamma_dT, dln_gamma_dn)

print *, "ln_gamma = ", ln_gamma
print *, "dln_gamma_dT = ", dln_gamma_dT
print *, "dln_gamma_dn = ", dln_gamma_dn
end program main
```

### Excess enthalpy \( (H^E) \)
From the Gibbs-Helmholtz equation [1]:

$$
 \left(\frac{\partial \left(\frac{G^E}{T} \right)}{\partial T} 
 \right)_P = \frac{-H^E}{T^2}
$$
We can calculate the excess enthalpy from the excess Gibbs energy as:

$$
 H^E = G^E - T \frac{\partial G^E}{\partial T}
$$

The derivatives of the excess enthalpy can be calculated as:

$$
\frac{\partial H^E}{\partial T} =
-T \frac{\partial^2 G^E}{\partial T^2}
$$

$$
\frac{\partial H^E}{\partial n_i} = \frac{\partial G^E}{\partial n_i}
- T \frac{\partial^2 G^E}{\partial T \partial n_i}
$$

#### Example
To simplify the example, we will use the same code as the previous example,
ommited here for brevity. The following code calculates the excess enthalpy and
its derivatives.

```fortran
real(pr) :: H_E, dH_E_dT, dH_E_dn(nc)

! Calculate excess enthalpy and its derivatives
call model%excess_enthalpy(n, T, HE, dHE_dT, dHE_dn)

print *, "HE = ", HE
print *, "dHE_dT = ", dHE_dT
print *, "dHE_dn = ", dHE_dn
```


### Excess entropy \( (S^E) \)
Finally, excess entropy can be calculated from the excess Gibbs energy and
excess enthalpy as:

$$
S^E = \frac{H^E - G^E}{T}
$$

The derivatives of the excess entropy can be calculated as:

$$
\frac{\partial S^E}{\partial T} = \frac{(\frac{\partial H^E}{\partial T} 
- \frac{\partial G^E}{\partial T})T - H^E + G^E}{T^2}
$$

$$
\frac{\partial S^E}{\partial n_i} = \frac{\frac{\partial H^E}
{\partial n_i} - \frac{\partial G^E}{\partial n_i}}{T}
$$

#### Example
To simplify the example, we will use the same code as the previous example,
ommited here for brevity. The following code calculates the excess entropy and
its derivatives.

```fortran
real(pr) :: SE, dSE_dT, dSE_dn(nc)

! Calculate excess entropy and its derivatives
call model%excess_entropy(n, T, SE, dSE_dT, dSE_dn)

print *, "S_E = ", SE
print *, "dSE_dT = ", dSE_dT
print *, "dSE_dn = ", dSE_dn
```


### References
[1] https://en.wikipedia.org/wiki/Gibbs%E2%80%93Helmholtz_equation