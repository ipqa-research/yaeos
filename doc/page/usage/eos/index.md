---
title: Equations of State
ordered_subpage: cubics
---

[TOC]

# Equations of State

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
thermodynamic properties can be calculated from the residual Helmholtz free
energy and its derivatives. `yaeos` is based on this approach.

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
A^r(\mathbf{n}, V, T) = -\int_{\infty}^{V} \left(P(\mathbf{n},V,T) -
\frac{nRT}{V} \right) dV
$$

We are ready to go now. We leave an example of how to instantiate a
Peng-Robinson EOS model in the next code block if you want to start evaluating
properties right away.

```fortran
program main
   use yaeos

   ! We will instantiate a Peng-Robinson EOS model with two components

   ! Set the variable `model` as a generic `ArModel`
   class(ArModel), allocatable :: model

   ! Set the variables that we're going to use as variable length arrays
   real(pr) :: tc(2), pc(2), w(2)

   tc = [190, 310]   ! Critical temperatures [K]
   pc = [14, 30]     ! Critical pressures [bar]
   w = [0.001, 0.03] ! Acentric factors [-]

   ! Now we set our model as the PengRobinson76 equation of state.
   model = PengRobinson76(tc, pc, w)

end program main
```

The examples provided for each property must be treated as a continuation of
the previous code block (before the end program statement).

# Thermodynamic properties

The reference used for the equations of each property is "Thermodynamic Models:
Fundamentals and Computational Aspects" by Michael L. Michelsen and Jørgen M.
Mollerup - 2th edition (abbreviated as: MM). For each property, we are going to
refer to the chapter and equation for a rapid search.

The MM book sometimes express the equations in terms of \(F\) (reduced
residual Helmholtz free energy):

$$
    F(\mathbf{n}, V, T) = \frac{A^r(\mathbf{n}, V, T)}{RT}
$$

On this documentation section we are going to express all the properties in
terms of \(A^r(\mathbf{n}, V, T)\) so a little algebraic work is required.

Finally, we are going to use the MM notation, for that:

$n$ is the total number of moles, and \(n_i\) is the number of moles of
component \(i\). Therefore:

$$
n = \sum_i n_i
$$

The moles number vector is denoted as \(\mathbf{n}\).

On partial derivatives, the subscript indicates the variables that are kept constant. For example:

$$
    \left(\frac{\partial P}{\partial V} \right)_{V,\mathbf{n}}
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
- \frac{nR}{V}
$$

$$
\left(\frac{\partial P}{\partial n_i} \right)_{V,T} =
-\left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_{T} + \frac{RT}{V}
$$

```fortran
pressure: block
    real(pr) :: n(2), T, V, P, dPdV, dPdT, dPdn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 0.1_pr           ! Set volume to 0.1 L

    ! Calculate pressure and its derivatives
    call model%pressure(n, V, T, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)

    print *, "Pressure: ", P
    print *, "dPdV: ", dPdV
    print *, "dPdT: ", dPdT
    print *, "dPdn: ", dPdn

end block pressure
```

### Volume: \(V(\mathbf{n}, P, T)\)
As you must know, given a temperature and pressure, the volume has not a unique
solution. In the case of a cubic EOS, there are three possible solutions for
the volume at a given temperature, pressure and composition. For this reason,
the volume is calculated iteratively. Provided \(\mathbf{n}\), \(P\), and \(T\)
we fix \(\mathbf{n}\) and \(T\) and iterate over the volume until the pressure
is the same as the input pressure.

To learm how the initial guess for \(V\) is obtained, please refer to the book
"Thermodynamic Models: Fundamentals and Computational Aspects" by Michael L.
Michelsen and Jørgen M. Mollerup. Volume derivatives can be calculated in terms
of pressure derivatives, this will appear later on this documentation.

```fortran
volume: block
    real(pr) :: n(2), T, P, V

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set pressure to 1 bar

    ! Calculate different volume roots
    call model%volume(n, P, T, V=V, root="liquid")
    print *, "Liquid volume: ", V

    call model%volume(n, P, T, V=V, root="vapor")
    print *, "Vapor volume: ", V

    call model%volume(n, P, T, V=V, root="stable")
    print *, "Stable volume root: ", V

end block volume
```

## Residual properties at \((\mathbf{n}, V, T)\)

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

There is a relation between both residual properties. In this section we are
going to calculate all the \((\mathbf{n},V,T)\) residual properties.

### Fugacity coefficients (V,T): \(\ln \phi_i (V,T)\)

**MM - Chapter 2 - Equations 13, 14, 15, and 16**

Natural logarithm of fugacity coefficients specifing \(V\) and \(T\) are
calculated as follows:

$$
\ln \hat{\phi}_i = \frac{1}{RT} \left( \frac{\partial A^r}{\partial n_i}
\right)_{V,T} - \ln Z
$$

Remember that the compressibility factor \(Z\) is calculated as:

$$
Z = \frac{PV}{nRT}
$$

Next, the derivatives:

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
 \frac{1}{RT} \frac{\left(\frac{\partial P}{\partial n_i} 
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
    real(pr) :: n(2), T, V, lnPhi(2), dlnPhidP(2), dlnPhidT(2), dlnPhidn(2,2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    call eos%lnphi_vt(&
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
    real(pr) :: n(2), T, V, lnPhi(2), P, dPdV, dPdT, dPdn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    call eos%lnphi_vt(&
      n, V, T, lnPhi=lnPhi, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
      )

end block lnphi_vt_p
```



### Fugacity (V,T): \(\ln f_i (V,T)\)

Alternative way of calculating fugacity directly from the residual Helmholtz
free energy:

```fortran
fugacity_vt: block
    real(pr) :: n(2), T, V, lnf(2), dlnfdV(2), dlnfdT(2), dlnfdn(2,2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    call eos%lnfug_vt(&
      n, V, T, lnf, &
      dlnfdV, dlnfdT, dlnfdn, &
      )

end block fugacity_vt
```

### Residual entropy (V,T): \(S^r(\mathbf{n},V,T)\)

**MM - Chapter 2 - Equation 17**

We start explaining the residual entropy because it is the easiest property to
calculate. Also, helps us to understand how to calculate the next properties.

The residual entropy is calculated as follows:

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
    real(pr) :: n(2), T, V, S, dSdV, dSdT, dSdn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    ! Calculate residual entropy and its derivatives
    call eos%residual_entropy(n, V, T, Sr=Sr, SrT=SrT, SrV=SrV, Srn=Srn)

    print *, "Residual entropy: ", Sr
    print *, "SrV: ", SrV
    print *, "SrT: ", SrT
    print *, "Srn: ", Srn

end block residual_entropy
```

### Residual enthalpy (V,T): \(H^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation XII)**

**MM - Chapter 2 - Equation 20**

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
    real(pr) :: n(2), T, V, Hr, HrV, HrT, Hrn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    ! Calculate residual enthalpy and its derivatives
    call eos%residual_enthalpy(n, V, T, Hr=Hr, HrV=HrV, HrT=HrT, Hrn=Hrn)

    print *, "Residual enthalpy: ", Hr
    print *, "HrV: ", HrV
    print *, "HrT: ", HrT
    print *, "Hrn: ", Hrn

end block residual_enthalpy
```

### Residual Gibbs free energy (V,T): \(G^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation VI)**

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
    real(pr) :: n(2), T, V, Gr, GrV, GrT, Grn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    ! Calculate residual Gibbs free energy and its derivatives
    call eos%residual_gibbs(n, V, T, Gr=Gr, GrV=GrV, GrT=GrT, Grn=Grn)

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
    real(pr) :: n(2), T, V, Ur, UrV, UrT, Urn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    ! Calculate residual internal energy and its derivatives
    call eos%internal_energyresidual(n, V, T, Ur=Ur, UrV=UrV, UrT=UrT, Urn=Urn)

    print *, "Residual internal energy: ", Ur
    print *, "UrV: ", UrV
    print *, "UrT: ", UrT
    print *, "Urn: ", Urn

end block residual_internal_energy
```

### Residual constant volume heat capacity (V,T): \(C_V^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation III)**

The residual constant volume heat capacity is calculated as follows:

$$
C_V^r(\mathbf{n},V,T) = - T \left(\frac{\partial^2 A^r}{\partial T^2}
\right)_{V,\mathbf{n}}
$$

```fortran
residual_cv: block
    real(pr) :: n(2),T, V, Cv

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    ! Calculate residual constant volume heat capacity
    call eos%residual_cv(n, V, T, Cv=Cv)

    print *, "Residual constant volume heat capacity: ", Cv

end block residual_cv
```

### Residual constant pressure heat capacity (V,T): \(C_P^r(\mathbf{n},V,T)\)

**MM - Chapter 1 - Table 6 (Equation VII)**

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
terms of residual Helmholtz (\mathbf{n},V,T).

$$
C_P^r(\mathbf{n},V,T) = - T \left(\frac{\partial^2 A^r}{\partial T^2}
\right)_{V,\mathbf{n}} - T \; \frac{\left(- \left(\frac{\partial^2
A^r}{\partial V \partial T} \right)_{V,\mathbf{n}} +
\frac{nR}{V}\right)^2}{\left(-\left(\frac{\partial^2 A^r}{\partial V^2}
\right)_{T,\mathbf{n}} - \frac{nRT}{V^2}\right)} - nR
$$


```fortran
residual_cp: block
    real(pr) :: n(2), T, V, Cp

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    V = 1.0_pr           ! Set volume to 1 liter

    ! Calculate residual constant pressure heat capacity
    call eos%residual_cp(n, V, T, Cp=Cp)

    print *, "Residual constant pressure heat capacity: ", Cp

end block residual_cp
```

## Residual properties at \((\mathbf{n}, P, T)\)

As explained in the "Residual properties at \((\mathbf{n}, V, T)\)" section, a residual
property is the difference between the actual property and the ideal gas
property. At \((\mathbf{n}, P, T)\), the difference is made against an ideal gas at
the same pressure of the real fluid, for that, the real fluid and the ideal
gas are at different volumes:

$$
    M^r(\mathbf{n},P,T) = M(\mathbf{n},P,T) - M^{ig}(\mathbf{n},P,T)
$$

Some properties change because of that, other stays the same. To obtain the
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
    V}\right)_{T,\mathbf{n}}} = \frac{-\left(\frac{\partial^2 A^r}{\partial V
    \partial n_i} \right)_T + \frac{RT}{V}}{-\left(\frac{\partial^2
    A^r}{\partial V^2} \right)_{T,\mathbf{n}} - \frac{nRT}{V^2}}
$$

From **MM - Chapter 1 - Equation A36**

$$
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} = -
    \frac{\left(\frac{\partial P}{\partial T}
    \right)_{V,\mathbf{n}}}{\left(\frac{\partial P}{\partial V}
    \right)_{T,\mathbf{n}}} = - \frac{- \left(\frac{\partial^2 A^r}{\partial V
    \partial T} \right)_{\mathbf{n}} + \frac{nR}{V}}{-\left(\frac{\partial^2
    A^r}{\partial V^2} \right)_{T,\mathbf{n}} - \frac{nRT}{V^2}}
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
    real(pr) :: n(2), T, P, lnPhi(2), dlnPhidP(2), dlnPhidT(2), dlnPhidn(2,2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%lnphi_pt(&
      n, P, T, root_type="liquid", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )

    call eos%lnphi_pt(&
      n, P, T, root_type="vapor", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )
    
    call eos%lnphi_pt(&
      n, P, T, root_type="stable", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
      )

end block lnphi_pt
```

As in the \(\ln \phi_i (V,T)\) method you can take advantage and the retrieve
the value of the calculated volume and the pressure derivatives.

```fortran
lnphi_pt_v: block
    real(pr) :: n(2), T, P, V, lnPhi(2), dPdV, dPdT, dPdn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%lnphi_pt(&
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
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} - nR \; ln \; Z -
    n R T \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    T}\right)_{P,\mathbf{n}} - \frac{1}{T}\right)
$$

$$
    \left(\frac{\partial A^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial A^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial A^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} - RT \; ln \; Z - nRT
    \left(\frac{1}{V} \left(\frac{\partial V}{\partial n_i}\right)_{P,T} -
    \frac{1}{n}\right)
$$

```fortran
ar_pt: block
    real(pr) :: n(2), T, P, Ar, ArP, ArT, Arn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%helmholtz_residual_pt(&
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
    \left(\frac{\partial S^r(\mathbf{n},P,T)}{\partial P}\right)_{T,n} =
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial P}\right)_{T,n} + nR
    \left(\frac{1}{P} + \frac{1}{V} \left(\frac{\partial V}{\partial
    P}\right)_{T,n} \right)
$$

$$
    \left(\frac{\partial S^r(\mathbf{n},P,T)}{\partial T}\right)_{P,n} =
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial T}\right)_{V,n} +
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial T}\right)_{P,n} + nR \left(\frac{1}{V}
    \left(\frac{\partial V}{\partial T}\right)_{P,n} - \frac{1}{T}\right)
$$

$$
    \left(\frac{\partial S^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial S^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} + R \; ln \; Z + nR
    \left(\frac{1}{V} \left(\frac{\partial V}{\partial n_i}\right)_{P,T} -
    \frac{1}{n}\right)
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
    \left(\frac{\partial H^r(\mathbf{n},P,T)}{\partial P}\right)_{T,n} =
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial P}{\partial V}\right)^{-1}_{T,n}
$$

$$
    \left(\frac{\partial H^r(\mathbf{n},P,T)}{\partial T}\right)_{P,n} =
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial T}\right)_{V,n} +
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial T}\right)_{P,n}
$$

$$
    \left(\frac{\partial H^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial H^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T}
$$

```fortran
hr_pt: block
    real(pr) :: n(2), T, P, Hr, HrP, HrT, Hrn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%enthalpy_residual_pt(&
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

Therefore, the derivatives are the same as the residual Helholtz free energy
\((\mathbf{n}, P, T)\):

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
    \left(\frac{\partial V}{\partial T}\right)_{P,\mathbf{n}} - nR \; ln \; Z -
    n R T \left(\frac{1}{V} \left(\frac{\partial V}{\partial
    T}\right)_{P,\mathbf{n}} - \frac{1}{T}\right)
$$

$$
    \left(\frac{\partial G^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial G^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial G^r(\mathbf{n},V,T)}{\partial V}\right)_{T,\mathbf{n}}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T} - RT \; ln \; Z - nRT
    \left(\frac{1}{V} \left(\frac{\partial V}{\partial n_i}\right)_{P,T} -
    \frac{1}{n}\right)
$$

```fortran
gr_pt: block
    real(pr) :: n(2), T, P, Gr, GrP, GrT, Grn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%gibbs_residual_pt(&
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

Therefore, its derivatives are the same as enthalpy \((\mathbf{n}, P, T)\):

$$
    \left(\frac{\partial U^r(\mathbf{n},P,T)}{\partial P}\right)_{T,n} =
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial P}{\partial V}\right)^{-1}_{T,n}
$$

$$
    \left(\frac{\partial U^r(\mathbf{n},P,T)}{\partial T}\right)_{P,n} =
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial T}\right)_{V,n} +
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial T}\right)_{P,n}
$$

$$
    \left(\frac{\partial U^r(\mathbf{n},P,T)}{\partial n_i}\right)_{P,T} =
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial n_i}\right)_{V,T} +
    \left(\frac{\partial U^r(\mathbf{n},V,T)}{\partial V}\right)_{T,n}
    \left(\frac{\partial V}{\partial n_i}\right)_{P,T}
$$

```fortran
ur_pt: block
    real(pr) :: n(2), T, P, Ur, UrP, UrT, Urn(2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%gibbs_residual_pt(&
       n, P, T, root_type="stable", Ur=Ur, UrP=UrP, UrT=UrT, Urn=Urn &
       )
end block ur_pt
```

### Residual constant volume heat capacity (P,T): \(C_V^r(\mathbf{n},P,T)\)

**MM - Chapter 1 - Table 6 (Equation X)**

Residual constant volume heat capacity \((\mathbf{n}, P, T)\) can be obtained
as follows:

$$
    C_v(\mathbf{n}, P, T) = C_v(\mathbf{n}, V, T)
$$

```fortran
cv_pt: block
    real(pr) :: n(2), T, P, Cv

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%Cv_residual_pt(n, P, T, root_type="stable", Cv=Cv)
end block cv_pt
```

### Residual constant pressure heat capacity (P,T): \(C_P^r(\mathbf{n},P,T)\)

**MM - Chapter 1 - Table 6 (Equation XIV)**

Residual constant pressure heat capacity \((\mathbf{n}, P, T)\) can be obtained
as follows:

$$
    C_p(\mathbf{n}, P, T) = C_p(\mathbf{n}, V, T)
$$

```fortran
cp_pt: block
    real(pr) :: n(2), T, P, Cp

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%Cp_residual_pt(n, P, T, root_type="stable", Cp=Cp)
end block cp_pt
```

## Excess properties

### Activity coefficients

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
    \left(\frac{\partial \ln \; \gamma_i}{\partial T}\right)_{P,n} =
    \left(\frac{\partial \ln \; \hat{\phi}_i}{\partial T}\right)_{P,n} -
    \left(\frac{\partial \ln \; \phi_i}{\partial T}\right)_{P}
$$

$$
    \left(\frac{\partial \ln \; \gamma_i}{\partial P}\right)_{T,n} =
    \frac{\overline{V}_i(\mathbf{n}, P, T) - v_i(P, T)}{RT}
$$

$$
    \left(\frac{\partial \ln \; \gamma_i}{\partial n_j}\right)_{P,T} =
    \left(\frac{\partial \ln \; \hat{\phi}_i}{\partial n_j}\right)_{P,T}
$$

Where:

$$
\overline{V}_i(\mathbf{n}, P, T) = \left(\frac{\partial V}{\partial
n_i}\right)_{P,T} = - \frac{\left(\frac{\partial P}{\partial n_i}\right)_{V,T}}{\left(\frac{\partial P}{\partial V}\right)_{T,\mathbf{n}}}
$$

and \(v_i(P, T)\) is the partial molar volume of component \(i\) in the pure
state at the same pressure and temperature of the mixture.

```fortran
ln_gamma: block
    real(pr) :: n(2), T, P, lngamma(2)
    real(pr) :: dlngammadT(2), dlngammadP(2), dlngammadn(2,2)

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%ln_activity_coefficient(&
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
    \left(\frac{\partial G^E}{\partial n_j}\right)_{P,T} = R T \sum_i n_i
    \left(\frac{\partial ln \; \gamma_i}{\partial n_j}\right)_{P,T} + R T
    \; ln \; \gamma_j
$$

```fortran
gibbs_excess: block
    real(pr) :: n(2), T, P, GE, dGEdT(2), dGEdP(2), dGE_dn(2,2)
    
    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%gibbs_excess(&
       n, P, T, root_type="stable", GE=GE, dGEdT=dGEdT, dGEdP=dGEdP, &
       dGE_dn=dGE_dn &
       )
end block gibbs_excess
```

### Excess Gibbs free energy (P, T): \(H^E(\mathbf{n}, P, T)\)

**MM - Chapter 1 - Table 9 (Equation III)**

The excess enthalpy can be calculated in terms of the activity coefficients as
follows:

$$
    H^E(\mathbf{n}, P, T) = -R T^2 \sum_i n_i \left(\frac{\partial ln \;
    \gamma_i}{\partial T}\right)_{P,\mathbf{n}}
$$

```fortran
enthalpy_excess: block
    real(pr) :: n(2), T, P, HE

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%enthalpy_excess(n, P, T, root_type="stable", HE=HE)
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
    real(pr) :: n(2), T, P, VE
    
    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%volume_excess(&
       n, P, T, root_type="stable", VE=VE &
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
    real(pr) :: n(2), T, P, SE

    n = [3.0_pr, 7.0_pr] ! Number of moles of each component [mol]
    T = 300.0_pr         ! Set temperature to 300 K
    P = 1.0_pr           ! Set volume to 1 bar

    call eos%entropy_excess(n, P, T, root_type="stable", SE=SE)
end block entropy_excess
```

### Excess Helmholtz (P, T): \(A^E(\mathbf{n}, P, T)\)

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

# References
[1] 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
Fundamentals & computational aspects (2. ed). Tie-Line Publications.
