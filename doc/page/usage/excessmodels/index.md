---
title: Gibbs Excess Models
ordered_subpage: nrtl.md
                 nrtl-hv.md
                 uniquac.md
                 unifac_models
---

[TOC]

# Gibbs Excess Models

## What problem does an excess-Gibbs model solve?

Mixing water with an alcohol does not simply dilute two otherwise unchanged
liquids: each molecule encounters new neighbors. An excess-Gibbs model describes
how these unlike contacts change the chemical potential (the driving force for
transfer or separation) relative to an ideal liquid solution. In `yaeos`, these
models extend [[GeModel]] and calculate activity coefficients and excess thermal
properties from temperature and amounts. They do **not** calculate liquid density
or vapor properties on their own.

An **excess** property compares a real mixture with an **ideal solution** at the
same temperature, pressure, and composition:

$$
    G^E = G - G^{I,S}
$$

$$
    H^E = H - H^{I,S} = \Delta H_{mix}
$$

$$
    S^E = S - S^{I,S}
$$

Unlike excess enthalpy, excess Gibbs energy and excess entropy are **not** the
total changes on mixing. An ideal solution already has an entropy of mixing:

$$
\Delta G_{mix}=RT\sum_i n_i\ln x_i+G^E,\qquad
\Delta S_{mix}=-R\sum_i n_i\ln x_i+S^E.
$$

A **residual** property in an equation of state instead measures departure from
an ideal gas reference. It is not interchangeable with an excess property. Also,
the word "residual" inside UNIQUAC/UNIFAC names just their energetic contribution
to the *excess* Gibbs energy, not an ideal-gas departure. The
[thermodynamics essentials page](../fundamentals.html) develops this distinction
in more depth, including why an excess quantity cannot be added to a residual
one.

## Activity, activity coefficient, and fugacity

With a pure-liquid standard state, the dimensionless activity is
\(a_i=x_i\gamma_i\), and
\(\mu_i=\mu_i^{*,L}+RT\ln a_i\). The activity coefficient \(\gamma_i\)
corrects ideal dilution: it is 1 in an ideal solution. A value above 1 raises
the chemical potential relative to that ideal solution; below 1 lowers it.
It is **not** a mole fraction and need not lie between 0 and 1.

Fugacity \(f_i\) is a pressure-like representation of chemical potential. In a
liquid, \(f_i^L=x_i\gamma_i f_i^{*,L}\); in a vapor described by an EOS,
\(f_i^V=y_i\varphi_i P\). Here \(\varphi_i\) is a **fugacity coefficient**, not
an activity coefficient. Vapor-liquid equilibrium requires equal fugacities,
so activity coefficients alone are not a complete vapor-liquid calculation:
pure-component reference fugacities and a vapor model are also needed. A
liquid-liquid calculation with the same standard states can cancel those
reference terms, but still requires stability and phase-split calculations.

## Choosing a model

| Available information / goal | Starting point |
|:--|:--|
| Reliable fitted binary liquid-mixture parameters | [NRTL](nrtl.html) or [UNIQUAC](uniquac.html) |
| Molecular structures, but no fitted molecular pair parameters | [UNIFAC family](unifac_models/index.html), provided all required groups and interactions are covered |
| An excess model intended for Huron-Vidal cubic-EOS coupling | [NRTL-HV](nrtl-hv.html), with matching EOS covolumes and fitted parameters |

Parameterized models can represent the data used to fit them well, but the
model name alone supplies no parameters. Group-contribution models reuse
interactions between chemical fragments across molecules; they are predictive
estimates, not a guarantee of accuracy for a new system. Choose using measured
data near your temperature/composition and the intended equilibrium type, not
just whether a constructor accepts the inputs.

## Inputs, units, and shared methods

- `n(:)` contains component **mole amounts**, with `n >= 0` and `sum(n) > 0`.
  Component order must match all model vectors and matrix indices. A mole-fraction
  vector is valid as a one-mole basis; mass fractions are not interchangeable.
- `T` is absolute temperature in **kelvin**, never degrees Celsius.
- `ln_activity_coefficient` returns **natural logarithms**; use `exp` to get
  \(\gamma\). Optional derivatives are with respect to `T` or mole amounts,
  not independently varied mole fractions.
- The gas constant [[yaeos__constants(module)]] uses bar L/(mol K). Thus extensive
  `Ge` and `He` are in **bar L**, and `Se` and `Cpe` are in **bar L/K** when `n`
  is in mol. Multiply by 100 to convert bar L to joules. Divide by `sum(n)` for
  molar properties.
- Scaling every amount by the same factor leaves activity coefficients unchanged
  and scales extensive properties by that factor. Mole-number derivatives depend
  on the chosen amount basis (the `dlngammadn` matrix scales inversely).

These implementations have no explicit pressure argument in their excess-property
methods. They cannot supply pressure-dependent liquid volumes or replace a full
EOS at high pressure. Each concrete model supplies `excess_gibbs`; the shared
methods derive the properties below from it.

Because the shared methods are inherited from [[GeModel]], the same calls work
for every model on this site. That uniformity is a convenience, not a statement
that the models are interchangeable: they differ in how they describe
nonideality, in how many parameters they need, and in where they are reliable.

## Worked example: reading a result

The complete program below evaluates an ethanol–water mixture with UNIFAC. It is
worth running once and reading the output carefully, because the numbers are
easily misread.

```fortran
program excess_properties
use yaeos, only: pr, Groups, setup_unifac, UNIFAC
implicit none

type(UNIFAC) :: model

integer, parameter :: nc = 2

type(Groups) :: molecules(nc)

real(pr) :: n(nc)
real(pr) :: T
real(pr) :: ln_gamma(nc)
real(pr) :: dln_gamma_dT(nc)
real(pr) :: dln_gamma_dn(nc, nc)
real(pr) :: Ge
real(pr) :: He
real(pr) :: Se
real(pr) :: Cpe

T = 298.15_pr
n = [0.5_pr, 0.5_pr]

! Ethanol, then water: output arrays follow this same order.
molecules(1)%groups_ids = [1, 2, 14]
molecules(1)%number_of_groups = [1, 1, 1]
molecules(2)%groups_ids = [16]
molecules(2)%number_of_groups = [1]

! setup UNIFAC model
model = setup_unifac(molecules)

! Calculate ln_gamma and its derivatives
call model%ln_activity_coefficient(n, T, ln_gamma, dln_gamma_dT, dln_gamma_dn)

print *, "ln_gamma = ", ln_gamma
print *, "dln_gamma_dT = ", dln_gamma_dT
print *, "dln_gamma_dn = ", dln_gamma_dn
print *, "gamma = ", exp(ln_gamma)
call model%excess_gibbs(n, T, Ge=Ge)
call model%excess_enthalpy(n, T, He=He)
call model%excess_entropy(n, T, Se=Se)
call model%excess_Cp(n, T, Cpe)
print *, "GE, HE [bar L] = ", Ge, He
print *, "SE, CpE [bar L/K] = ", Se, Cpe
print *, "GE - (HE - T*SE), near zero = ", Ge - (He - T*Se)
end program excess_properties
```

This is a complete program. The first output entry refers to ethanol and the
second to water. Negative `ln_gamma` is allowed; `exp(ln_gamma)` must be positive.
The final energy identity is a useful consistency check, but does not validate
the model against experiments. A positive `He` means mixing absorbs heat relative
to ideal-solution enthalpy at that temperature; `Cpe` is only the *excess*, not
the total heat capacity.

To interpret the output, work through it in this order.

1. **`ln_gamma` is a logarithm.** `exp(ln_gamma)` gives the coefficient. A
   coefficient above 1 means the component is *less* stable in the mixture than
   ideal dilution would predict; below 1 means more stable. Neither the
   logarithm nor the coefficient is a mole fraction, and neither is bounded by 1.
2. **Multiply by the mole fraction for an activity.** `x_i * gamma_i` is the
   dimensionless activity, the quantity that appears in the chemical potential.
   A large `gamma_i` for a component with a small `x_i` may still mean a small
   activity.
3. **The energy outputs are extensive.** `Ge` and `He` are in bar·L for the
   whole `n` vector (1 mol here, so numerically molar); `Se` and `Cpe` are in
   bar·L/K. Multiply by 100 for J and J/K. Divide by `sum(n)` for molar values.
   A positive `Ge` means mixing is less favorable than ideal mixing; a positive
   `He` means mixing absorbs heat.
4. **`dln_gamma_dn` is a mole-number derivative**, so it scales inversely with
   the amount basis: multiplying `n` by ten divides these entries by ten, while
   `ln_gamma` itself is unchanged. `dln_gamma_dT` is unaffected by scaling.
5. **The identity check.** `Ge - (He - T*Se)` must vanish at the level of
   round-off, because \(G^E=H^E-TS^E\). A nonzero value means a bug in the
   model or in your call, not a physical finding.

What this example does *not* tell you is whether the prediction is accurate.
That requires comparing with measured data for the same mixture.

## Routines and Methods
### \( \ln \gamma_i \)

Activity coefficients are calculated from the excess Gibbs energy using the
following expression:

$$
    \ln \gamma_i = \left(\frac{\partial \left(\frac{G^E}{RT} \right)}{\partial n_i}\right)_{T,P,n_{j\ne i}}
$$

In words: the activity coefficient is the sensitivity of the excess Gibbs energy
to adding a little more of component `i`. A mixture whose excess Gibbs energy
responds strongly to extra `i` has a large `ln gamma_i`. Because `G^E` is
extensive and the derivative is taken at fixed amounts of the other
components, `ln gamma_i` depends only on composition and temperature, not on how
much material is present.

#### Example

The complete program in [Worked example](#worked-example-reading-a-result) shows
the call in context. The essential part is:

```fortran
call model%ln_activity_coefficient(n, T, ln_gamma, dln_gamma_dT, dln_gamma_dn)
```

Pass `n` (mole amounts) and `T` (kelvin). The three outputs are optional; asking
only for `lngamma` is cheaper because the derivative kernels are skipped.

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


The complete program that produces these numbers is in
[Worked example](#worked-example-reading-a-result).

### Excess enthalpy \( (H^E) \)
From the Gibbs-Helmholtz equation [1]:

$$
 \left(\frac{\partial \left(\frac{G^E}{T} \right)}{\partial T} 
 \right)_{P,\mathbf{n}} = \frac{-H^E}{T^2}
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
The following is an **extension snippet**, not a separate program. Put declarations
with the declarations of the complete example, and calls before its `end program`.

```fortran
real(pr) :: HE_value
real(pr) :: dHE_dT
real(pr) :: dHE_dn(nc)

! Calculate excess enthalpy and its derivatives
call model%excess_enthalpy(n, T, HE_value, dHE_dT, dHE_dn)

print *, "HE = ", HE_value
print *, "dHE_dT = ", dHE_dT
print *, "dHE_dn = ", dHE_dn
```


### Excess entropy \( (S^E) \)
Finally, excess entropy can be calculated from the excess Gibbs energy and
excess enthalpy as:

$$
S^E = \frac{H^E - G^E}{T}
$$

Replacing \(H^E\):

$$
S^E= \frac{G^E - T \frac{\partial G^E}{\partial T} - G^E}{T} = 
-\frac{\partial G^E}{\partial T}
$$

The derivatives of the excess entropy can be calculated as:

$$
\frac{\partial S^E}{\partial T} = -\frac{\partial^2 G^E}{\partial T^2}
$$

$$
\frac{\partial S^E}{\partial n_i} = 
-\frac{\partial^2 G^E}{\partial T \partial n_i}
$$

#### Example
This is another **extension snippet**; separate its declarations and executable
statements when inserting it in the complete example.

```fortran
real(pr) :: SE_value
real(pr) :: dSE_dT
real(pr) :: dSE_dn(nc)

! Calculate excess entropy and its derivatives
call model%excess_entropy(n, T, SE_value, dSE_dT, dSE_dn)

print *, "S_E = ", SE_value
print *, "dSE_dT = ", dSE_dT
print *, "dSE_dn = ", dSE_dn
```

### Excess heat capacity \( (C_{p}^{E}) \)

Excess heat capacity can be calculated as:

$$
C_p^E = \left(\frac{\partial H^E}{\partial T}\right)_{\mathbf{n}} = -T
\frac{\partial^2 G^E}{\partial T^2}
$$

#### Example

Extension snippet (declaration belongs above the first executable statement):

```fortran
real(pr) :: C_p_E

! Calculate excess heat capacity
call model%excess_Cp(n, T, C_p_E)

print *, "C_p_E = ", C_p_E
```

## Checks, limits, and next steps

Before using a phase-equilibrium result, work through the inputs once more.
Confirm the component order matches the model, that the temperature is inside
the range the parameters were fitted for, that the signs and units of the
parameters are what you think they are, and, for the group-contribution models,
that every interaction your mixture needs is actually covered. Scale `n` and
confirm the activity coefficients do not change. In the pure-component limit the
coefficient of the **present** component should approach 1 for a consistently
referenced model; a trace solute need not have a coefficient of 1. Never
evaluate an all-zero amount vector.

A large activity coefficient does not by itself prove that a liquid splits.
Whether a feed is unstable is a separate question, answered by a stability
analysis and a phase-equilibrium calculation; see the
[phase-equilibrium guide](../phase_equilibrium/index.html). And a converged
prediction is not an accurate one: compare with measured data before drawing
conclusions from it.

For EOS coupling, continue to the
[cubic mixing rules](../eos/cubics/mixing.html). Standard NRTL, UNIQUAC, and
UNIFAC are not general electrolyte, reacting-mixture, or near-critical models,
and they supply no volumetric properties on their own.

### References
[1] https://en.wikipedia.org/wiki/Gibbs%E2%80%93Helmholtz_equation
