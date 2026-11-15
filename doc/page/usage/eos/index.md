---
title: Equations of State
ordered_subpage: cubics
                 gerg2008.md
                 pc_saft.md
---

[TOC]

# Equations of State

## Start here

An EOS answers questions such as “what pressure does this amount of fluid exert
in this volume?” and “how far is this fluid from ideal-gas behavior?”. These
are **single homogeneous-phase** calculations; evaluating an EOS at an overall
mixture composition does not automatically split it into equilibrium phases.
See [thermodynamic fundamentals](../fundamentals.html) for the broader workflow.

**For practical use:** read the units below, run the complete example, then
look up the property you need in the table of contents. The detailed equations
are an optional reference, not a prerequisite for calling the routines.

| What you know / need | Routine or section |
| --- | --- |
| Amounts, total volume, temperature → pressure | `pressure` |
| Amounts, pressure, temperature → volume of a chosen root | `volume` |
| Nonideality for phase-equilibrium calculations | `lnphi_pt` or `lnphi_vt` |
| Energy or heat-capacity correction to an ideal gas | `*_residual_pt` or `*_residual_vt` |
| Deviation from ideal-solution mixing | Excess properties |

### Units and amount basis

The default unit system uses mol, L, K, and bar, with
\(R = 0.08314462618\ \mathrm{bar\,L\,mol^{-1}\,K^{-1}}\).
Do not pass Pa or m³ without converting them.

- The Fortran array `n(:)` contains component **amounts**, \(\mathbf{n}=(n_1,\ldots,n_c)\),
  in mol; the scalar total is \(n=\sum_i n_i\).
- Mole fractions are \(z_i=n_i/n\). They sum to one and specify composition,
  not the amount of material. Passing `n=z` selects a one-mole basis.
- `V` is **total** volume in L; molar volume is \(v=V/n\) in L/mol.
  If you multiply every amount by ten, also multiply `V` by ten to keep
  the same density and pressure.
- `Ar`, `Hr`, `Gr`, and `Ur` are total energies in bar L, not molar energies;
  divide by `sum(n)` for bar L/mol. One bar L is 100 J.
  Entropy and heat capacities are in bar L/K. Derivative units follow from
  division by the differentiated variable, e.g. `dPdV` in bar/L.
- `lnPhi` and `lngamma` are dimensionless vectors, one entry per component.
  Composition-derivative matrices have entry `(i,j)` equal to the derivative
  of output `i` with respect to amount `n(j)`, not mole fraction `z(j)`.

Choose a model and parameters appropriate to your fluid and range; a successful
calculation is not evidence that the model is accurate there. See
[cubic models](cubics/index.html), [GERG-2008](gerg2008.html), and
[PC-SAFT](pc_saft.html) for model-specific inputs and limitations.

## Why Helmholtz energy?

`yaeos` is a library based on Equation-of-State calculations. An Equation of
State (EOS) is a mathematical model that describes the relationship between
thermodynamic properties of a system, typically pressure (P), volume (V), and
temperature (T). These equations are fundamental in physics, chemistry, and
engineering, especially in thermodynamics and fluid mechanics.

EOS models are used to predict phase behavior, thermodynamic properties, and
equilibrium conditions of gases, liquids, solids, and mixtures. They are widely
applied in areas such as chemical engineering, petroleum engineering, and
material science. From the ideal gas law, to cubic equations of state and SAFT,
to more complex models, EOS are essential tools for the study and design of
processes.

The classic way of defining an EOS is though a mathematical equation that
express the pressure as a function of the volume and temperature. For example,
the Van der Waals EOS is defined as:

$$
P = \frac{RT}{v - b} - \frac{a}{v^2}
$$

Where \(P\) is the pressure, \(v\) is the molar volume, \(T\) is the
temperature, \(R\) is the gas constant, and \(a\) and \(b\) are parameters that
depend on the substance. For mixtures, the \(a\) and \(b\) parameters are
calculated from the pure components parameters using a mixing rule.

The modern way of defining an EOS is through a mathematical model that express
the residual Helmholtz free energy as a function of the mole number, volume and
temperature. This approach is more general and allows for a more flexible
definition of the EOS in the context of implementing it in code. All
residual thermodynamic properties can be calculated from the residual Helmholtz
free energy and its derivatives. Total caloric properties also require an
ideal-gas contribution and a consistent reference state. `yaeos` is based on
this approach.

This section of the documentation will guide you through the thermodynamic
properties that can be calculated with yaeos using Helmholtz free energy
models. Additionally, this page serves as a summary of how all these properties
can be derived from the Helmholtz free energy and its derivatives.

First, you might be used to the classic \(P(\mathbf{n},V,T)\) way of defining an
EOS. This is normal since it is the most common way of explaining them in
undergraduate thermodynamics courses. So, a very reasonable question is: "how
can I obtain the expression for the residual Helmholtz free energy
\(A^r(\mathbf{n}, V, T)\) from the \(P(\mathbf{n}, V, T)\) expression?". The answer
is the following:

$$
 A^r(\mathbf{n}, V, T) = -\int_{\infty}^{V} \left(P(\mathbf{n},V',T) -
\frac{nRT}{V'} \right) dV'
$$

The integral holds at fixed amounts and temperature, with the residual tending
to zero in the dilute ideal-gas limit.

### A complete first calculation

This standalone program constructs a Peng–Robinson model and evaluates a vapor
root. The critical constants below are illustrative inputs, not a validated
parameter set for named substances.

```fortran
program main
   use yaeos, only: pr, R, ArModel, PengRobinson76
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   ! We will instantiate a Peng-Robinson EOS model with two components

   ! Set the variable `model` as a generic `ArModel`
   class(ArModel), allocatable :: model

   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: n(2)
   real(pr) :: T
   real(pr) :: P
   real(pr) :: V
   real(pr) :: Pcheck
   real(pr) :: Z
   real(pr) :: lnPhi(2)
   real(pr) :: Hr

   tc = [190.0_pr, 310.0_pr] ! Critical temperatures [K]
   pc = [14.0_pr, 30.0_pr]   ! Critical pressures [bar]
   w = [0.001_pr, 0.03_pr]   ! Acentric factors [-]

   ! Now we set our model as the PengRobinson76 equation of state.
   model = PengRobinson76(tc, pc, w)

   n = [0.3_pr, 0.7_pr] ! One mole, composition z = [0.3, 0.7]
   T = 300.0_pr
   P = 1.0_pr
   call model%volume(n, P, T, V=V, root_type="vapor")
   if (.not. ieee_is_finite(V)) error stop "Nonfinite volume"
   if (V <= 0.0_pr) error stop "Volume solve failed"
   call model%pressure(n, V, T, P=Pcheck)
   if (.not. ieee_is_finite(Pcheck)) error stop "Nonfinite pressure"
   if (abs(Pcheck - P) > 1.0e-6_pr) error stop "Pressure mismatch"
   call model%lnphi_vt(n, V, T, lnPhi=lnPhi)
   call model%enthalpy_residual_vt(n, V, T, Hr=Hr)
   Z = Pcheck*V/(sum(n)*R*T)
   print *, "Total volume [L]: ", V
   print *, "Molar volume [L/mol]: ", V/sum(n)
   print *, "Compressibility factor [-]: ", Z
   print *, "ln(phi_i) [-]: ", lnPhi
   print *, "Residual molar enthalpy [J/mol]: ", 100.0_pr*Hr/sum(n)

end program main
```

Here `Z=1` is ideal-gas behavior; deviations show how the pressure–volume
relation differs from an ideal gas. `lnPhi=0` means a unit fugacity coefficient,
not zero fugacity. A negative `Hr` means an enthalpy below the ideal-gas value
at the same temperature and composition; it is not an absolute enthalpy.

All subsequent code blocks are **snippets**: insert them after model
construction and before `end program main`, not as separate programs. Their
`block` constructs give local variables their own scope. Optional derivative
outputs may be omitted when only the property itself is needed.

# Thermodynamic properties

The reference used for the equations of each property is "Thermodynamic Models:
Fundamentals and Computational Aspects" by Michael L. Michelsen and Jørgen M.
Mollerup - 2nd edition (abbreviated as: MM). For each property, we are going to
refer to the chapter and equation for a rapid search.

The MM book sometimes express the equations in terms of \(F\) (reduced
residual Helmholtz free energy):

$$
    F(\mathbf{n}, V, T) = \frac{A^r(\mathbf{n}, V, T)}{RT}
$$

On this documentation section we are going to express all the properties in
terms of \(A^r(\mathbf{n}, V, T)\) so a little algebraic work is required.

Finally, we are going to use the MM notation, for that:

\(n\) is the total number of moles, and \(n_i\) is the number of moles of
component \(i\). Therefore:

$$
n = \sum_i n_i
$$

The moles number vector is denoted as \(\mathbf{n}\).

On partial derivatives, the subscript indicates the variables that are kept constant. For example:

$$
    \left(\frac{\partial P}{\partial V} \right)_{T,\mathbf{n}}
$$

means the partial derivative of pressure with respect to volume keeping
temperature and moles vector constant. With a partial derivative respect to the
mole number of component \(i\) we have:

$$
    \left(\frac{\partial P}{\partial n_i} \right)_{V,T}
$$

With a rigorous notation, the partial derivative of pressure with respect to
the mole number of component \(i\) keeping volume and temperature constant
should be:

$$
    \left(\frac{\partial P}{\partial n_i} \right)_{V,T} =
    \left(\frac{\partial P}{\partial n_i} \right)_{V,T,n_{j \neq i}}
$$

On the notation the subscript \(n_{j \neq i}\) means that all the mole numbers
of the other components are kept constant. This is omitted for simplicity.

## Pressure and volume

### Pressure: \(P(\mathbf{n}, V, T)\)

**MM - Chapter 2 - Equations 7, 9, 10, and 11**

Pressure can be calculated from the residual Helmholtz free energy as follows:

$$
P = - \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}} 
+ \frac{nRT}{V}
$$

$$
\left(\frac{\partial P}{\partial V} \right)_{T,\mathbf{n}} =
-\left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,\mathbf{n}} -
\frac{nRT}{V^2}
$$

$$
\left(\frac{\partial P}{\partial T} \right)_{V,\mathbf{n}} =
- \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_{\mathbf{n}} +
\frac{nR}{V}
$$

$$
\left(\frac{\partial P}{\partial n_i} \right)_{V,T} =
-\left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_{T} + \frac{RT}{V}
$$

```fortran
pressure: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: P
    real(pr) :: dPdV
    real(pr) :: dPdT
    real(pr) :: dPdn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate pressure and its derivatives
    call model%pressure(n, V, T, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)

    print *, "Pressure: ", P
    print *, "dPdV: ", dPdV
    print *, "dPdT: ", dPdT
    print *, "dPdn: ", dPdn

end block pressure
```

### Volume: \(V(\mathbf{n}, P, T)\)
At fixed amounts, pressure, and temperature, an EOS may have more than one
volume root. A conventional cubic EOS has one or three real algebraic roots
(counting multiplicity); not every root is physically admissible. Three
distinct physical phases are not implied by three roots. Other model families
need not follow this root count. `yaeos` uses numerical volume solvers; this is
an implementation choice, not a consequence of nonuniqueness alone.

`root_type="liquid"` seeks the dense root and `"vapor"` the dilute root.
They may return the same volume when only one root exists. `"stable"` selects
the lower-Gibbs-energy candidate at the **specified composition**; it does not
test stability against separation into phases of different compositions.
For equilibrium phase amounts and compositions, use the
[phase-equilibrium routines](../phase_equilibrium/index.html).

Check that a returned volume is finite and positive, and re-evaluate pressure
to check the solve, as in the first example. Root solvers can fail, especially
near limiting states. A mechanically stable homogeneous branch normally has
\((\partial P/\partial V)_{T,\mathbf{n}}<0\); this is necessary, not sufficient,
for full mixture stability.

To learn more about volume initialization, please refer to the book
"Thermodynamic Models: Fundamentals and Computational Aspects" by Michael L.
Michelsen and Jørgen M. Mollerup. Volume derivatives can be calculated in terms
of pressure derivatives, this will appear later on this documentation.

```fortran
volume: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: V

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    ! Calculate different volume roots
    call model%volume(n, P, T, V=V, root_type="liquid")
    print *, "Liquid volume: ", V

    call model%volume(n, P, T, V=V, root_type="vapor")
    print *, "Vapor volume: ", V

    call model%volume(n, P, T, V=V, root_type="stable")
    print *, "Stable volume root: ", V

end block volume
```

## Residual properties at (n, V, T)

As you may know, residual properties are the difference between the actual
property and the ideal gas property at the same temperature and volume (in this
case). For that, for a generic residual property \(M^r\), we have:

$$
    M^r(\mathbf{n},V,T) = M(\mathbf{n},V,T) - M^{\text{ig}}(\mathbf{n},V,T)
$$

Notice that the ideal gas property \(M^{\text{ig}}(\mathbf{n},V,T)\) is
calculated at the same volume of the real fluid. If the ideal gas is at the
same temperature and volume of the real fluid, the ideal gas is not at the same
pressure as the real fluid. For that, two different residual properties can be
calculated:

- \((\mathbf{n},V,T)\)
- \((\mathbf{n},P,T)\)

Here we call the same-\(V,T\) difference a **residual** and the same-\(P,T\)
difference a **departure**; the API uses `residual` in both routine names.
The sections below retain that API terminology. The difference is the ideal-gas
reference, not a different real-fluid state: when comparing the two, \(P\) is
the EOS pressure at \((\mathbf{n},V,T)\).

These energy, entropy, and heat-capacity corrections are extensive: scaling
\(\mathbf{n}\) and \(V\) together by a factor scales them by that factor.
Fugacity coefficients are intensive and are **not** differences between total
properties. Their conventional ideal-gas reference is at the same \(P,T\),
even when the routine takes \(V,T\) as inputs.

### Fugacity coefficients (V,T): \(\ln \phi_i (V,T)\)

**MM - Chapter 2 - Equations 13, 14, 15, and 16**

Fugacity is an effective pressure for chemical equilibrium. For a phase of
composition \(z_i=n_i/n\), \(f_i=z_i\hat\phi_i P\); \(\hat\phi_i\) corrects
the ideal-gas partial pressure. Equilibrium requires equal component fugacities
between phases, not equal fugacity coefficients. Logarithms make products
additive and help numerical solvers handle large fugacity ratios. For an ideal
gas, \(\ln\hat\phi_i=0\); positive and negative values mean coefficients
above and below one, respectively.

Natural logarithms of fugacity coefficients specifying \(V\) and \(T\) are
calculated as follows:

$$
\ln \hat{\phi}_i = \frac{1}{RT} \left( \frac{\partial A^r}{\partial n_i}
\right)_{V,T} - \ln Z
$$

Remember that the compressibility factor \(Z\) is calculated as:

$$
Z = \frac{PV}{nRT}
$$

The derivative outputs are useful for Newton solvers, stability calculations,
and property sensitivities. **Despite the `vt` input specification**,
`dlnPhidT`, `dlnPhidP`, and `dlnPhidn` are derivatives at the constraints shown
below: temperature at fixed \(P,\mathbf{n}\), pressure at fixed
\(T,\mathbf{n}\), and amounts at fixed \(P,T\). Pressure derivatives returned
by the same routine retain the constraints in the pressure section.

$$
\left(\frac{\partial \ln \hat{\phi_i}}{\partial T} \right)_{P,\mathbf{n}} = 
\frac{1}{RT} \left(\frac{\partial^2 A^r}{\partial T \partial n_i}\right)_V
- \frac{1}{RT^2} \left(\frac{\partial A^r}{\partial n_i}\right)_{V,T}
+ \frac{1}{T} + \frac{1}{RT} \frac{\left(\frac{\partial P}{\partial n_i} 
\right)_{V,T} \left(\frac{\partial P}{\partial T} \right)_{V,\mathbf{n}}}{
\left(\frac{\partial P}{\partial V}\right)_{T,\mathbf{n}}}
$$

$$
\left(\frac{\partial \ln \hat{\phi_i}}{\partial P} \right)_{T,\mathbf{n}} = 
 -\frac{1}{RT} \frac{\left(\frac{\partial P}{\partial n_i}
\right)_{V,T}}{\left(\frac{\partial P}{\partial V}\right)_{T,\mathbf{n}}} 
- \frac{1}{P}
$$

$$
\left(\frac{\partial \ln \hat{\phi_i}}{\partial n_j} \right)_{P,T} =
\frac{1}{n} + \frac{1}{RT}
\left(\frac{\partial^2 A^r}{\partial n_i \partial n_j}
+ \frac{\left(\frac{\partial P}{\partial n_i} \right)_{V,T}
\left(\frac{\partial P}{\partial n_j} \right)_{V,T}}
{\left(\frac{\partial P}{\partial V} \right)_{T,\mathbf{n}}}\right)
$$

```fortran
lnphi_vt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: lnPhi(2)
    real(pr) :: dlnPhidP(2)
    real(pr) :: dlnPhidT(2)
    real(pr) :: dlnPhidn(2,2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    call model%lnphi_vt(&
      n, V, T, lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )

end block lnphi_vt

```

The subroutine `lnphi_vt` also allows to obtain the value of pressure and its
derivatives at the specified volume and temperature. Since those values are
needed to calculate the fugacity coefficients, you can take advantage of this
feature to avoid calculating them twice.


```fortran

lnphi_vt_p: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: lnPhi(2)
    real(pr) :: P
    real(pr) :: dPdV
    real(pr) :: dPdT
    real(pr) :: dPdn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    call model%lnphi_vt(&
      n, V, T, lnPhi=lnPhi, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
      )

end block lnphi_vt_p
```



### Fugacity (V,T): \(\ln f_i (V,T)\)

`lnfug_vt` returns \(\ln(f_i/1\,\mathrm{bar})\), using fugacity numerically
in bar, rather than its coefficient. For positive component amounts,

$$
\ln(f_i/1\,\mathrm{bar}) = \ln(z_i) + \ln\hat\phi_i + \ln(P/1\,\mathrm{bar}).
$$

Unlike `lnphi_vt`, its derivatives are at fixed \(V\) or \(T\): `dlnfdT`
holds \(V,\mathbf{n}\) fixed, `dlnfdV` holds \(T,\mathbf{n}\) fixed, and
`dlnfdn` holds \(V,T\) and the other component amounts fixed. In particular,

$$
\left(\frac{\partial\ln f_i}{\partial V}\right)_{T,\mathbf{n}}
=-\frac{1}{RT}\left(\frac{\partial P}{\partial n_i}\right)_{V,T},
\qquad
\left(\frac{\partial\ln f_i}{\partial n_j}\right)_{V,T}
=\frac{1}{RT}\frac{\partial^2 A^r}{\partial n_i\partial n_j}
+\frac{\delta_{ij}}{n_i}.
$$

Use positive amounts in this example. At zero amount, \(\ln f_i\) has a
singular limit: the current routine sets the corresponding `lnf` entry to
zero as an implementation convention, **not** the physical log fugacity,
and composition derivatives can divide by zero.

```fortran
fugacity_vt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: lnf(2)
    real(pr) :: dlnfdV(2)
    real(pr) :: dlnfdT(2)
    real(pr) :: dlnfdn(2,2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    call model%lnfug_vt(&
      n, V, T, lnf=lnf, &
      dlnfdV=dlnfdV, dlnfdT=dlnfdT, dlnfdn=dlnfdn &
      )

end block fugacity_vt
```

### Residual entropy (V,T): \(S^r(\mathbf{n},V,T)\)

**MM - Chapter 2 - Equation 17**

We start explaining the residual entropy because it is the easiest property to
calculate. Also, helps us to understand how to calculate the next properties.

Entropy describes the thermal and mixing contribution to the free energy.
`Sr` is the total entropy correction relative to an ideal gas at the same
amounts, volume, and temperature (bar L/K). It is not total entropy, and its
sign is not a stability test. The residual entropy is calculated as follows:

$$
S^r = - \left(\frac{\partial A^r}{\partial T} \right)_{V,\mathbf{n}}
$$

And its derivatives:

$$
\left(\frac{\partial S^r}{\partial T} \right)_{V,\mathbf{n}} =
- \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,\mathbf{n}}
$$


$$
\left(\frac{\partial S^r}{\partial V} \right)_{T,\mathbf{n}} =
- \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_{\mathbf{n}}
$$

$$
\left(\frac{\partial S^r}{\partial n_i} \right)_{V,T} =
- \left(\frac{\partial^2 A^r}{\partial n_i \partial T} \right)_V
$$

```fortran
residual_entropy: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: Sr
    real(pr) :: SrV
    real(pr) :: SrT
    real(pr) :: Srn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate residual entropy and its derivatives
    call model%entropy_residual_vt(n, V, T, Sr=Sr, SrT=SrT, SrV=SrV, Srn=Srn)

    print *, "Residual entropy: ", Sr
    print *, "SrV: ", SrV
    print *, "SrT: ", SrT
    print *, "Srn: ", Srn

end block residual_entropy
```

### Residual enthalpy (V,T): \(H^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation XII)**

**MM - Chapter 2 - Equation 20**

Enthalpy \(H=U+PV\) combines internal energy and the pressure–volume term,
making it useful in constant-pressure energy balances. `Hr` is only its
nonideal correction, in bar L. To obtain total enthalpy, add a compatible
ideal-gas value, \(H=H^{ig}+H^r\); for a mixture,
\(H^{ig}=\sum_i n_i h_i^{ig}(T)\), with a specified reference for each species.
Because ideal-gas enthalpy depends only on temperature and amounts, the
same-\(V,T\) residual and same-\(P,T\) departure enthalpies coincide at the
same real-fluid state. Their derivatives have different constraints.

The residual enthalpy is calculated as follows:

$$
H^r(\mathbf{n},V,T) = H^r(\mathbf{n},P,T) = A^r(\mathbf{n},V,T) + T
S^r(\mathbf{n},V,T) + PV - nRT
$$

We know that pressure can be calculated as:

$$
P = - \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}} +
\frac{nRT}{V}
$$

Then, we can obtain:

$$
P - \frac{nRT}{V} = - \left(\frac{\partial A^r}{\partial V}
\right)_{T,\mathbf{n}}
$$

$$
PV - nRT = - V \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

And we also know how to calculate residual entropy as:

$$
S^r = - \left(\frac{\partial A^r}{\partial T} \right)_{V,\mathbf{n}}
$$

Then, we can obtain the residual enthalpy as:

$$
H^r = A^r(\mathbf{n},V,T) - T \left(\frac{\partial A^r}{\partial T}
\right)_{V,\mathbf{n}}
- V \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

And its derivatives:

$$
\left(\frac{\partial H^r}{\partial T} \right)_{V,\mathbf{n}} =
- T \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,\mathbf{n}}
- V \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_{\mathbf{n}}
$$

$$
\left(\frac{\partial H^r}{\partial V} \right)_{T,\mathbf{n}} =
- T \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_{\mathbf{n}}
- V \left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,\mathbf{n}}
$$

$$
\left(\frac{\partial H^r}{\partial n_i} \right)_{V,T} =
\left(\frac{\partial A^r}{\partial n_i} \right)_{V,T}
- T \left(\frac{\partial^2 A^r}{\partial T \partial n_i} \right)_V
- V \left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_T
$$

```fortran
residual_enthalpy: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: Hr
    real(pr) :: HrV
    real(pr) :: HrT
    real(pr) :: Hrn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate residual enthalpy and its derivatives
    call model%enthalpy_residual_vt(n, V, T, Hr=Hr, HrV=HrV, HrT=HrT, Hrn=Hrn)

    print *, "Residual enthalpy: ", Hr
    print *, "HrV: ", HrV
    print *, "HrT: ", HrT
    print *, "Hrn: ", Hrn

end block residual_enthalpy
```

### Residual Gibbs free energy (V,T): \(G^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation VI)**

Gibbs energy \(G=H-TS\) is the potential used for equilibrium at fixed
pressure and temperature. The routine here returns a **same-volume** residual,
not the same-pressure departure used in root comparisons. In particular,
\(RT\sum_i n_i\ln\hat\phi_i\) is the same-pressure departure Gibbs energy;
subtract \(nRT\ln Z\) from the result below to obtain it.

The residual Gibbs free energy is calculated as follows:

$$
G^r(\mathbf{n},V,T) = A^r(\mathbf{n},V,T) + PV - nRT
$$

As with residual enthalpy we can easily deduce:

$$
PV - nRT = - V \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

For that:

$$
G^r(\mathbf{n},V,T) = A^r(\mathbf{n},V,T) - V \left(\frac{\partial
A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

And its derivatives:

$$
\left(\frac{\partial G^r}{\partial T} \right)_{V,\mathbf{n}} =
\left(\frac{\partial A^r}{\partial T} \right)_{V,\mathbf{n}}
- V \left(\frac{\partial^2 A^r}{\partial T \partial V} \right)_{\mathbf{n}}
$$

$$
\left(\frac{\partial G^r}{\partial V} \right)_{T,\mathbf{n}} =
- V \left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,\mathbf{n}}
$$

$$
\left(\frac{\partial G^r}{\partial n_i} \right)_{V,T} =
\left(\frac{\partial A^r}{\partial n_i} \right)_{V,T}
- V \left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_{T}
$$

```fortran
residual_gibbs: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: Gr
    real(pr) :: GrV
    real(pr) :: GrT
    real(pr) :: Grn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate residual Gibbs free energy and its derivatives
    call model%gibbs_residual_vt(n, V, T, Gr=Gr, GrV=GrV, GrT=GrT, Grn=Grn)

    print *, "Residual Gibbs free energy: ", Gr
    print *, "GrV: ", GrV
    print *, "GrT: ", GrT
    print *, "Grn: ", Grn

end block residual_gibbs
```

Little check of algebraic transformations:

From **MM - Chapter 2 - Equation 22** we know that:

$$
    G^r(\mathbf{n},P,T) = H^r(\mathbf{n},P,T) - T S^r(\mathbf{n},P,T)
$$

From **MM - Chapter 1 - Table 6 (Equations IX, XII, and XIII)** we know that:

$$
    S^r(\mathbf{n},P,T) = S^r(\mathbf{n},V,T) + nR \; ln Z
$$

$$
    H^r(\mathbf{n},P,T) = H^r(\mathbf{n},V,T)
$$

$$
    G^r(\mathbf{n},P,T) = G^r(\mathbf{n},V,T) - nRT \; ln Z
$$

Replacing on the first equation we have:

$$
G^r(\mathbf{n},V,T) - nRT \; ln Z = H^r(\mathbf{n},V,T) - T
\left(S^r(\mathbf{n},V,T) + nR \; ln Z \right)
$$

Which of course lead us to:

$$
G^r(\mathbf{n},V,T) = H^r(\mathbf{n},V,T) - T S^r(\mathbf{n},V,T)
$$

And we have established that:

$$
G^r(\mathbf{n},V,T) = A^r(\mathbf{n},V,T) - V \left(\frac{\partial
A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

$$
S^r = - \left(\frac{\partial A^r}{\partial T} \right)_{V,\mathbf{n}}
$$

$$
H^r = A^r(\mathbf{n},V,T) - T \left(\frac{\partial A^r}{\partial T}
\right)_{V,\mathbf{n}}
- V \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

Replacing all the properties in the equation lead us to:

$$ 
A^r(\mathbf{n},V,T) - V \left(\frac{\partial A^r}{\partial V}
\right)_{T,\mathbf{n}} = A^r(\mathbf{n},V,T) - T \left(\frac{\partial
A^r}{\partial T} \right)_{V,\mathbf{n}}
- V \left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}} - T \left(
- \left(\frac{\partial A^r}{\partial T} \right)_{V,\mathbf{n}}\right)
$$

With a very little algebra the equality is satisfied:

$$ 
A^r(\mathbf{n},V,T) - V \left(\frac{\partial A^r}{\partial V}
\right)_{T,\mathbf{n}} = A^r(\mathbf{n},V,T) - V \left(\frac{\partial
A^r}{\partial V} \right)_{T,\mathbf{n}}
$$

### Residual internal energy (V,T): \(U^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation IV)**

Internal energy \(U\) excludes the \(PV\) term in enthalpy. `Ur` is its
total nonideal correction (bar L), to be added to ideal-gas internal energy.
The residual internal energy is calculated as follows:

$$
U^r(\mathbf{n},V,T) = A^r(\mathbf{n},V,T) + T S^r(\mathbf{n},V,T)
$$

Therefore:

$$
U^r = A^r - T \left(\frac{\partial A^r}{\partial T} \right)_{V,\mathbf{n}}
$$

And its derivatives:

$$
\left(\frac{\partial U^r}{\partial T} \right)_{V,\mathbf{n}} =
- T \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,\mathbf{n}}
$$

$$
\left(\frac{\partial U^r}{\partial V} \right)_{T,\mathbf{n}} =
\left(\frac{\partial A^r}{\partial V} \right)_{T,\mathbf{n}}
- T \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_{\mathbf{n}}
$$

$$
\left(\frac{\partial U^r}{\partial n_i} \right)_{V,T} =
\left(\frac{\partial A^r}{\partial n_i} \right)_{V,T}
- T \left(\frac{\partial^2 A^r}{\partial T \partial n_i} \right)_{V}
$$

```fortran
residual_internal_energy: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: Ur
    real(pr) :: UrV
    real(pr) :: UrT
    real(pr) :: Urn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate residual internal energy and its derivatives
    call model%internal_energy_residual_vt(n, V, T, Ur=Ur, UrV=UrV, UrT=UrT, Urn=Urn)

    print *, "Residual internal energy: ", Ur
    print *, "UrV: ", UrV
    print *, "UrT: ", UrT
    print *, "Urn: ", Urn

end block residual_internal_energy
```

### Residual constant volume heat capacity (V,T): \(C_V^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation III)**

Heat capacity measures the energy change per kelvin on a specified path.
\(C_V=(\partial U/\partial T)_{V,\mathbf{n}}\); the `Cv` output here is
\(C_V^r=C_V-C_V^{ig}\), not total \(C_V\), in bar L/K.
The residual constant volume heat capacity is calculated as follows:

$$
C_V^r(\mathbf{n},V,T) = - T \left(\frac{\partial^2 A^r}{\partial T^2}
\right)_{V,\mathbf{n}}
$$

```fortran
residual_cv: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: Cv

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate residual constant volume heat capacity
    call model%Cv_residual_vt(n, V, T, Cv=Cv)

    print *, "Residual constant volume heat capacity: ", Cv

end block residual_cv
```

### Residual constant pressure heat capacity (V,T): \(C_P^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation VII)**

\(C_P=(\partial H/\partial T)_{P,\mathbf{n}}\), even when its state is
specified by \(V,T\). The `Cp` output is \(C_P^r=C_P-C_P^{ig}\), not total
\(C_P\). Add \(\sum_i n_i c_{P,i}^{ig}(T)\) in the same units for the total;
a negative residual heat capacity alone does not mean a negative total heat
capacity. In particular, `HrT` from `enthalpy_residual_vt` is **not**
\(C_P^r\), since it holds volume rather than pressure fixed.

The residual constant pressure heat capacity is calculated as follows:

$$
C_P^r(\mathbf{n},V,T) = C_V^r(\mathbf{n},V,T) + T \left(\frac{\partial
P}{\partial T}\right)_{V,\mathbf{n}} \left(\frac{\partial V}{\partial
T}\right)_{P,\mathbf{n}} - nR
$$

From **MM - Chapter 1 - Equation A37** we know that:

$$ 
\left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} = -
\frac{\left(\frac{\partial P}{\partial T}
\right)_{V,\mathbf{n}}}{\left(\frac{\partial P}{\partial V}
\right)_{T,\mathbf{n}}}
$$

With a direct replacement:

$$
C_P^r(\mathbf{n},V,T) = C_V^r(\mathbf{n},V,T) - T \frac{\left(\frac{\partial P}
{\partial T} \right)^2_{V,\mathbf{n}}}{\left(\frac{\partial P}
{\partial V} \right)_{T,\mathbf{n}}} - nR
$$

We can replace all the terms to obtain an expression completely expressed in
terms of residual Helmholtz \(A^r(\mathbf{n},V,T)\).

$$
C_P^r(\mathbf{n},V,T) = - T \left(\frac{\partial^2 A^r}{\partial T^2}
\right)_{V,\mathbf{n}} - T \; \frac{\left(- \left(\frac{\partial^2
A^r}{\partial V \partial T} \right)_{\mathbf{n}} +
\frac{nR}{V}\right)^2}{\left(-\left(\frac{\partial^2 A^r}{\partial
V^2} \right)_{T,\mathbf{n}} - \frac{nRT}{V^2}\right)} - nR
$$


```fortran
residual_cp: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: V
    real(pr) :: Cp

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 25.0_pr          ! Set total volume to 25 L (molar volume 25 L/mol)

    ! Calculate residual constant pressure heat capacity
    call model%Cp_residual_vt(n, V, T, Cp=Cp)

    print *, "Residual constant pressure heat capacity: ", Cp

end block residual_cp
```

The constant-pressure formulas divide by \((\partial P/\partial V)_{T,\mathbf{n}}\).
Near a spinodal or critical point this derivative approaches zero, so heat
capacities and transformed derivatives may become very large or numerically
ill-conditioned. They describe a homogeneous branch, not the heat capacity of
a two-phase sample whose phase amounts change with temperature.

## Residual properties at (n, P, T)

As explained in the "Residual properties at \((\mathbf{n}, V, T)\)" section, a residual
property is the difference between the actual property and the ideal gas
property. At \((\mathbf{n}, P, T)\), the difference is made against an ideal gas at
the same pressure of the real fluid, for that, the real fluid and the ideal
gas are at different volumes:

$$
    M^r(\mathbf{n},P,T) = M(\mathbf{n},P,T) - M^{ig}(\mathbf{n},P,T)
$$

For the same real-fluid state, \(H^r\), \(U^r\), \(C_V^r\), and \(C_P^r\)
are unchanged because their ideal-gas counterparts depend only on temperature
and amounts. \(A^r\), \(G^r\), and \(S^r\) require the \(\ln Z\) corrections
given below. For these three, a `_pt` routine is therefore not simply the
`_vt` residual evaluated at the solved volume.

To obtain the
expression of these properties in term of the \(\mathbf{n},V,T\) properties we will need
the derivatives of \(ln \; Z\):

$$
    ln \; Z = ln \; P + ln \; V - \ln (nRT)
$$

$$
    \left(\frac{\partial \, ln \; Z}{\partial P}\right)_{T,\mathbf{n}} =
    \frac{1}{P} + \frac{1}{V} \left(\frac{\partial V}{\partial
    P}\right)_{T,\mathbf{n}}
$$

$$
    \left(\frac{\partial \, ln \; Z}{\partial T}\right)_{P,\mathbf{n}} =
    \frac{1}{V} \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} -
    \frac{1}{T}
$$

$$
    \left(\frac{\partial \, ln \; Z}{\partial n_i}\right)_{P,T} =
    \frac{1}{V} \left(\frac{\partial V}{\partial n_i}\right)_{P,T} 
    - \frac{1}{n}
$$

If desired, all the volume derivatives can be expressed in terms of
\(A^r(\mathbf{n}, V, T)\) as:

$$
    \left(\frac{\partial V}{\partial P}\right)_{T,\mathbf{n}} =
    \frac{1}{\left(\frac{\partial P}{\partial V}\right)_{T,\mathbf{n}}} =
    \frac{1}{-\left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,\mathbf{n}}
    - \frac{nRT}{V^2}}
$$

From **MM - Chapter 1 - Equation A39**

$$
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} =
    -\frac{\left(\frac{\partial P}{\partial
    n_i}\right)_{V,T}}{\left(\frac{\partial P}{\partial
    V}\right)_{T,\mathbf{n}}} = -\frac{-\left(\frac{\partial^2 A^r}{\partial V
    \partial n_i} \right)_T + \frac{RT}{V}}{-\left(\frac{\partial^2
    A^r}{\partial V^2} \right)_{T,\mathbf{n}} - \frac{nRT}{V^2}}
$$

From **MM - Chapter 1 - Equation A36**

$$
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} = -
    \frac{\left(\frac{\partial P}{\partial T}
    \right)_{V,\mathbf{n}}}{\left(\frac{\partial P}{\partial V}
    \right)_{T,\mathbf{n}}} = - \frac{- \left(\frac{\partial^2 A^r}{\partial V
    \partial T} \right)_{\mathbf{n}} +
    \frac{nR}{V}}{-\left(\frac{\partial^2 A^r}{\partial V^2}
    \right)_{T,\mathbf{n}} - \frac{nRT}{V^2}}
$$

Finally, we have to understand the chain rule...

$$
    \left(\frac{\partial A^r(\mathbf{n}, V, T)}{\partial
    T}\right)_{P,\mathbf{n}} = \left(\frac{\partial A^r(\mathbf{n}, V,
    T)}{\partial T}\right)_{V,\mathbf{n}} + \left(\frac{\partial
    A^r(\mathbf{n}, V, T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}}
$$


### Fugacity coefficients (P,T): \(\ln \phi_i (P,T)\)

There is no need for a direct way of calculating natural logarithm the fugacity
coefficients specifying \(P\) and \(T\) since we can solve volume from those
specifications and then calculate the fugacity coefficients using the \(\ln
\phi_i (V,T)\) method. For that reason you can choose the root of the volume
that you want to use to calculate the fugacity coefficients in the same way as
the volume method.

```fortran

lnphi_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: lnPhi(2)
    real(pr) :: dlnPhidP(2)
    real(pr) :: dlnPhidT(2)
    real(pr) :: dlnPhidn(2,2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%lnphi_pt(&
      n, P, T, root_type="liquid", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )

    call model%lnphi_pt(&
      n, P, T, root_type="vapor", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )
    
    call model%lnphi_pt(&
      n, P, T, root_type="stable", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )

end block lnphi_pt
```

As in the \(\ln \phi_i (V,T)\) method you can take advantage and the retrieve
the value of the calculated volume and the pressure derivatives.

```fortran
lnphi_pt_v: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: V
    real(pr) :: lnPhi(2)
    real(pr) :: dPdV
    real(pr) :: dPdT
    real(pr) :: dPdn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%lnphi_pt(&
      n, P, T, V=V, root_type="stable", lnPhi=lnPhi, &
      dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
      )

end block lnphi_pt_v

```

### Residual Helmholtz free energy (P, T): \(A^r(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 6 (Equation VIII)**

The \(A^r(\mathbf{n}, P, T)\) can be calculated as follows:

$$
    A^r(\mathbf{n}, P, T) = A^r(\mathbf{n}, V, T) - nRT \; ln \; Z
$$

And its derivatives:

$$
    \left(\frac{\partial \, A^r(\mathbf{n},P,T)}{\partial
    P}\right)_{T,\mathbf{n}} = \left(\frac{\partial \,
    A^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}} \left(\frac{\partial
    V}{\partial P}\right)_{T,\mathbf{n}} - n R T \left(\frac{1}{P} +
    \frac{1}{V} \left(\frac{\partial V}{\partial P}\right)_{T,\mathbf{n}}
    \right)
$$

$$
    \left(\frac{\partial \, A^r(\mathbf{n},P,T)}{\partial
    T}\right)_{P,\mathbf{n}} = \left(\frac{\partial A^r(\mathbf{n}, V,
    T)}{\partial T}\right)_{V,\mathbf{n}} + \left(\frac{\partial
    A^r(\mathbf{n}, V, T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} - nR \;
    ln \; Z - n R T \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    T}\right)_{P,\mathbf{n}} - \frac{1}{T}\right)
$$

$$
    \left(\frac{\partial A^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial A^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial A^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} - RT \; ln \; Z -
    nRT \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    n_i}\right)_{P,T} - \frac{1}{n}\right)
$$

```fortran
ar_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Ar
    real(pr) :: ArP
    real(pr) :: ArT
    real(pr) :: Arn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%helmholtz_residual_pt(&
       n, P, T, root_type="stable", Ar=Ar, ArP=ArP, ArT=ArT, Arn=Arn &
       )
end block ar_pt

```

### Residual entropy (P, T): \(S^r(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 6 (Equation IX)**

The residual entropy at \(\mathbf{n}, P, T\) can be calculated as follows:

$$
    S^r(\mathbf{n}, P, T) = S^r(\mathbf{n}, V, T) + nR \; \ln \; Z
$$

And its derivatives:

$$
    \left(\frac{\partial S^r(\mathbf{n},P,T)}{\partial P}\right)_{T,\mathbf{n}}
    = \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial
    V}\right)_{T,\mathbf{n}} \left(\frac{\partial V}{\partial
    P}\right)_{T,\mathbf{n}} + nR \left(\frac{1}{P} + \frac{1}{V}
    \left(\frac{\partial V}{\partial P}\right)_{T,\mathbf{n}} \right)
$$

$$
    \left(\frac{\partial S^r(\mathbf{n},P,T)}{\partial T}\right)_{P,\mathbf{n}}
    = \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial
    T}\right)_{V,\mathbf{n}} + \left(\frac{\partial
    S^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}} \left(\frac{\partial
    V}{\partial T}\right)_{P,\mathbf{n}} + nR \left(\frac{1}{V}
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} -
    \frac{1}{T}\right)
$$

$$
    \left(\frac{\partial S^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} + R \; ln \; Z +
    nR \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    n_i}\right)_{P,T} - \frac{1}{n}\right)
$$


### Residual Enthalpy (P, T): \(H^r(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 6 (Equation XII)**

In this case we have a property that doesn't change its value with the
specification.

$$
    H^r(\mathbf{n},P,T) = H^r(\mathbf{n},V,T)
$$

For that, its derivatives can be obtained as:

$$
    \left(\frac{\partial H^r(\mathbf{n},P,T)}{\partial P}\right)_{T,\mathbf{n}}
    = \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial
    V}\right)_{T,\mathbf{n}} \left(\frac{\partial P}{\partial
    V}\right)^{-1}_{T,\mathbf{n}}
$$

$$
    \left(\frac{\partial H^r(\mathbf{n},P,T)}{\partial T}\right)_{P,\mathbf{n}}
    = \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial
    T}\right)_{V,\mathbf{n}} + \left(\frac{\partial
    H^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}} \left(\frac{\partial
    V}{\partial T}\right)_{P,\mathbf{n}}
$$

$$
    \left(\frac{\partial H^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T}
$$

```fortran
hr_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Hr
    real(pr) :: HrP
    real(pr) :: HrT
    real(pr) :: Hrn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%enthalpy_residual_pt(&
       n, P, T, root_type="stable", Hr=Hr, HrP=HrP, HrT=HrT, Hrn=Hrn &
       )
end block hr_pt
```

### Residual Gibbs free energy (P, T): \(G^r(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 6 (Equation XIII)**

Residual Gibbs free energy \((\mathbf{n}, P, T)\) can be obtained as follows:

$$
    G^r(\mathbf{n}, P, T) = G^r(\mathbf{n}, V, T) - n R T \; \ln \; Z
$$

Therefore, the derivatives follow the same chain-rule pattern as residual
Helmholtz free energy \((\mathbf{n}, P, T)\), with \(G^r\) replacing \(A^r\):

$$
    \left(\frac{\partial \, G^r(\mathbf{n},P,T)}{\partial
    P}\right)_{T,\mathbf{n}} = \left(\frac{\partial \,
    G^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}} \left(\frac{\partial
    V}{\partial P}\right)_{T,\mathbf{n}} - n R T \left(\frac{1}{P} +
    \frac{1}{V} \left(\frac{\partial V}{\partial P}\right)_{T,\mathbf{n}}
    \right)
$$

$$
    \left(\frac{\partial \, G^r(\mathbf{n},P,T)}{\partial
    T}\right)_{P,\mathbf{n}} = \left(\frac{\partial G^r(\mathbf{n}, V,
    T)}{\partial T}\right)_{V,\mathbf{n}} + \left(\frac{\partial
    G^r(\mathbf{n}, V, T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} - nR \;
    ln \; Z - n R T \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    T}\right)_{P,\mathbf{n}} - \frac{1}{T}\right)
$$

$$
    \left(\frac{\partial G^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial G^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial G^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} - RT \; ln \; Z -
    nRT \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    n_i}\right)_{P,T} - \frac{1}{n}\right)
$$

```fortran
gr_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Gr
    real(pr) :: GrP
    real(pr) :: GrT
    real(pr) :: Grn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%gibbs_residual_pt(&
       n, P, T, root_type="stable", Gr=Gr, GrP=GrP, GrT=GrT, Grn=Grn &
       )
end block gr_pt
```

### Residual internal energy (P, T): \(U^r(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 6 (Equation XI)**

Residual internal energy \((\mathbf{n}, P, T)\) can be obtained as follows:

$$
    U^r(\mathbf{n}, P, T) = U^r(\mathbf{n}, V, T)
$$

Therefore, its derivatives follow the same chain-rule pattern as enthalpy
\((\mathbf{n}, P, T)\), with \(U^r\) replacing \(H^r\):

$$
    \left(\frac{\partial U^r(\mathbf{n},P,T)}{\partial P}\right)_{T,\mathbf{n}}
    = \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial
    V}\right)_{T,\mathbf{n}} \left(\frac{\partial P}{\partial
    V}\right)^{-1}_{T,\mathbf{n}}
$$

$$
    \left(\frac{\partial U^r(\mathbf{n},P,T)}{\partial T}\right)_{P,\mathbf{n}}
    = \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial
    T}\right)_{V,\mathbf{n}} + \left(\frac{\partial
    U^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}} \left(\frac{\partial
    V}{\partial T}\right)_{P,\mathbf{n}}
$$

$$
    \left(\frac{\partial U^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T}
$$

```fortran
ur_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Ur
    real(pr) :: UrP
    real(pr) :: UrT
    real(pr) :: Urn(2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%internal_energy_residual_pt(&
       n, P, T, root_type="stable", Ur=Ur, UrP=UrP, UrT=UrT, Urn=Urn &
       )
end block ur_pt
```

### Residual constant volume heat capacity (P,T): \(C_V^r(\mathbf{n},P,T)\)

**MM - Chapter 1 - Table 6 (Equation X)**

Residual constant volume heat capacity \((\mathbf{n}, P, T)\) can be obtained
as follows:

$$
    C_V^r(\mathbf{n}, P, T) = C_V^r(\mathbf{n}, V, T)
$$

```fortran
cv_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Cv

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%Cv_residual_pt(n, P, T, root_type="stable", Cv=Cv)
end block cv_pt
```

### Residual constant pressure heat capacity (P,T): \(C_P^r(\mathbf{n},P,T)\)

**MM - Chapter 1 - Table 6 (Equation XIV)**

Residual constant pressure heat capacity \((\mathbf{n}, P, T)\) can be obtained
as follows:

$$
    C_P^r(\mathbf{n}, P, T) = C_P^r(\mathbf{n}, V, T)
$$

```fortran
cp_pt: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Cp

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%Cp_residual_pt(n, P, T, root_type="stable", Cp=Cp)
end block cp_pt
```

## Excess properties (n, P, T)

Excess properties measure departure from an **ideal solution**, not from an
ideal gas. They are not generally equal to the full change on mixing: an ideal
solution already has entropy and Gibbs energy of mixing. At the same \(P,T\)
and composition,

$$
    M^E(\mathbf{n}, P, T) = M(\mathbf{n}, P, T) - M^{ideal}(\mathbf{n}, P, T)
$$

For \(M=H,U,V\), the ideal-solution property is the sum of pure-component
properties:

$$
    M^{ideal}(\mathbf{n}, P, T) = \sum_i n_i \; m_i(P, T)
$$

Where \(m_i(P, T)\) is the molar property of the pure component \(i\) at the
same pressure and temperature as the mixture.

For Gibbs and Helmholtz energies, add \(RT\sum_i n_i\ln z_i\) to this pure
sum; for entropy, add \(-R\sum_i n_i\ln z_i\). Those ideal-mixing terms
must be included before interpreting \(G^E\), \(A^E\), or \(S^E\).
For example, positive \(V^E\) means expansion relative to the pure-volume sum,
and negative \(H^E\) means an enthalpy below the pure-enthalpy sum.

The current routines use `root_type` for the **mixture** but always use
`"stable"` roots for the pure-component references. This matters if you want
a liquid-standard-state activity coefficient while a pure component's stable
state at the specified \(P,T\) is vapor. These routines do not automatically
construct hypothetical pure-liquid reference states. Derivatives assume a
smooth branch and are not derivatives across a reference-root switch.

### Activity coefficients (P, T): \(ln \gamma_i(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 9 (Equation 1)**

The natural logarithm of the activity coefficient is given by:

$$
    ln \; \gamma_i(\mathbf{n}, P, T) = ln \; \hat{\phi}_i(\mathbf{n}, P, T) -
    ln \; \phi_i(P, T)
$$

Where:

- \(\hat{\phi}_i(\mathbf{n}, P, T)\) is the fugacity coefficient of component
  \(i\) in the mixture

- \(\phi_i(P, T)\) is the fugacity coefficient of component \(i\) in the pure
  state at the same pressure and temperature of the mixture.

And its derivatives:

**MM - Chapter 1 - Table 9 (Equations V and VI)**

$$
    \left(\frac{\partial \ln \; \gamma_i}{\partial T}\right)_{P,\mathbf{n}} =
    \left(\frac{\partial \ln \; \hat{\phi}_i}{\partial T}\right)_{P,\mathbf{n}}
    - \left(\frac{\partial \ln \; \phi_i}{\partial T}\right)_{P}
$$

$$
    \left(\frac{\partial \ln \; \gamma_i}{\partial P}\right)_{T,\mathbf{n}} =
    \frac{\overline{V}_i(\mathbf{n}, P, T) - v_i(P, T)}{RT}
$$

$$
    \left(\frac{\partial \ln \; \gamma_i}{\partial n_j}\right)_{P,T} =
    \left(\frac{\partial \ln \; \hat{\phi}_i}{\partial n_j}\right)_{P,T}
$$

Where:

$$
    \overline{V}_i(\mathbf{n}, P, T) = \left(\frac{\partial V}{\partial
    n_i}\right)_{P,T} = - \frac{\left(\frac{\partial P}{\partial
    n_i}\right)_{V,T}}{\left(\frac{\partial P}{\partial
    V}\right)_{T,\mathbf{n}}}
$$

and \(v_i(P, T)\) is the molar volume of pure component \(i\) at the same
pressure and temperature as the mixture.

```fortran
ln_gamma: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: lngamma(2)
    real(pr) :: dlngammadT(2)
    real(pr) :: dlngammadP(2)
    real(pr) :: dlngammadn(2,2)

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%ln_activity_coefficient(&
       n, P, T, root_type="stable", &
       lngamma=lngamma, dlngammadT=dlngammadT, dlngammadP=dlngammadP, &
       dlngammadn=dlngammadn &
       )
end block ln_gamma
```

### Excess Gibbs free energy (P, T): \(G^E(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 9 (Equation II)**

The excess Gibbs free energy can be calculated in terms of the activity
coefficients as follows:

$$
    G^E(\mathbf{n}, P, T) = R T \sum_i n_i \; ln \; \gamma_i
$$

and its derivatives:

$$
    \left(\frac{\partial G^E}{\partial T}\right)_{P,\mathbf{n}} = R \sum_i n_i
    \; ln \; \gamma_i + R T \sum_i n_i \left(\frac{\partial ln \;
    \gamma_i}{\partial T}\right)_{P,\mathbf{n}}
$$

$$
    \left(\frac{\partial G^E}{\partial P}\right)_{T,\mathbf{n}} = R T \sum_i n_i
    \left(\frac{\partial ln \; \gamma_i}{\partial P}\right)_{T,\mathbf{n}}
$$

$$
    \left(\frac{\partial G^E}{\partial n_i}\right)_{P,T} = R T \; ln \; \gamma_i
$$

`Ge`, `GeP`, and `GeT` are scalars; `Gen` has one entry per component.
In particular, `Gen(i)` is a first amount derivative (bar L/mol), not a
matrix of second derivatives. `GeP` has units L and `GeT` bar L/K.

```fortran
gibbs_excess: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Ge
    real(pr) :: GeP
    real(pr) :: GeT
    real(pr) :: Gen(2)
    
    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%gibbs_excess(&
       n, P, T, root_type="stable", Ge=Ge, GeP=GeP, GeT=GeT, Gen=Gen &
       )
end block gibbs_excess
```

### Excess Enthalpy (P, T): \(H^E(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 9 (Equation III)**

The excess enthalpy can be calculated in terms of the activity coefficients as
follows:

$$
    H^E(\mathbf{n}, P, T) = -R T^2 \sum_i n_i \left(\frac{\partial ln \;
    \gamma_i}{\partial T}\right)_{P,\mathbf{n}}
$$

```fortran
enthalpy_excess: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: He

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%enthalpy_excess(n, P, T, root_type="stable", He=He)
end block enthalpy_excess
```

### Excess volume (P, T): \(V^E(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 9 (Equation IV)**

The excess volume can be calculated in terms of the activity coefficients as
follows:

$$
    V^E(\mathbf{n}, P, T) = R T \sum_i n_i \left(\frac{\partial ln \;
    \gamma_i}{\partial P}\right)_{T,\mathbf{n}}
$$

```fortran
volume_excess: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Ve
    
    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%volume_excess(&
       n, P, T, root_type="stable", Ve=Ve &
       )
end block volume_excess
```

### Excess entropy (P, T): \(S^E(\mathbf{n}, P, T)\)

The excess entropy can be calculated in terms of the activity coefficients as
follows:

$$
    S^E(\mathbf{n}, P, T) = \frac{H^E - G^E}{T}
$$

replacing \(H^E\) and \(G^E\)

$$
     S^E(\mathbf{n}, P, T) = \frac{-R T^2 \sum_i n_i \left(\frac{\partial ln \;
    \gamma_i}{\partial T}\right)_{P,\mathbf{n}} - R T \sum_i n_i \; ln \;
    \gamma_i}{T}
$$

$$
     S^E(\mathbf{n}, P, T) = -R T \sum_i n_i \left(\frac{\partial ln \;
    \gamma_i}{\partial T}\right)_{P,\mathbf{n}} - R \sum_i n_i \; ln \;
    \gamma_i
$$

```fortran
entropy_excess: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Se

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%entropy_excess(n, P, T, root_type="stable", Se=Se)
end block entropy_excess
```

### Excess Helmholtz free energy (P, T): \(A^E(\mathbf{n}, P, T)\)

**MM - Chapter 5 - Equation 1**

The excess Helmholtz free energy can be calculated in terms of the activity
coefficients as follows:

$$
    A^E(\mathbf{n}, P, T) = G^E(\mathbf{n}, P, T) - P V^E
$$

$$
    A^E(\mathbf{n}, P, T) = R T \sum_i n_i \; ln \; \gamma_i - P R T \sum_i n_i
    \left(\frac{\partial ln \; \gamma_i}{\partial P}\right)_{T,\mathbf{n}}
$$

```fortran
helmholtz_excess: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Ae

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%helmholtz_excess(n, P, T, root_type="stable", Ae=Ae)
end block helmholtz_excess
```

### Excess Internal Energy (P, T): \(U^E(\mathbf{n}, P, T)\)

**MM - Chapter 5 - Equation 4**

The excess internal energy can be calculated in terms of the activity
coefficients as follows:

$$
    U^E(\mathbf{n}, P, T) = A^E(\mathbf{n}, P, T) + T S^E(\mathbf{n}, P, T)
$$

Replacing:

$$
    U^E(\mathbf{n}, P, T) = R T \sum_i n_i \; ln \; \gamma_i - P R T \sum_i n_i
    \left(\frac{\partial ln \; \gamma_i}{\partial P}\right)_{T,\mathbf{n}} + T
    \left(-R T \sum_i n_i \left(\frac{\partial ln \; \gamma_i}{\partial
    T}\right)_{P,\mathbf{n}} - R \sum_i n_i \; ln \; \gamma_i \right)
$$

And with very little algebra:

$$
    U^E(\mathbf{n}, P, T) = - P R T \sum_i n_i
    \left(\frac{\partial ln \; \gamma_i}{\partial P}\right)_{T,\mathbf{n}}
    - R T^2 \sum_i n_i \left(\frac{\partial ln \; \gamma_i}{\partial
    T}\right)_{P,\mathbf{n}}
$$

```fortran
internal_energy_excess: block
    real(pr) :: n(2)
    real(pr) :: T
    real(pr) :: P
    real(pr) :: Ue

    n = [0.3_pr, 0.7_pr] ! Component amounts [mol]; total 1 mol, z = [0.3, 0.7]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    call model%internal_energy_excess(n, P, T, root_type="stable", Ue=Ue)
end block internal_energy_excess
```

# References
[1] Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
Fundamentals & computational aspects (2. ed). Tie-Line Publications.
