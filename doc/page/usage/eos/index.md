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

First of all, you might be used to the classic \(P(n,V,T)\) way of defining an
EOS. This is normal since it is the most common way of explaining them in
undergraduate thermodynamics courses. So, a very reasonable question is: "how
can I obtain the expression for the residual Helmholtz free energy
\(A^r(T,V,n)\) from the \(P(n,V,T)\) expression?". The answer is the following:

$$
A^r(n, V, T) = -\int_{\infty}^{V} \left(P(n,V,T) - \frac{nRT}{V} \right) dV
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
   real(pr), allocatable :: n(:), tc(:), pc(:), w(:)

   n = [0.3, 0.7]    ! Number of moles of each component [mol]
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
## Pressure: \(P(n, V, T)\)

Pressure can be calculated from the residual Helmholtz free energy as follows:

$$
P = - \left(\frac{\partial A^r}{\partial V} \right)_{T,n} 
+ \frac{nRT}{V}
$$

$$
\left(\frac{\partial P}{\partial V} \right)_{T,n} = 
-\left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,n} - \frac{nRT}{V^2}
$$

$$
\left(\frac{\partial P}{\partial T} \right)_{V,n} =
- \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_n + \frac{nR}{V}
$$

$$
\left(\frac{\partial P}{\partial n_i} \right)_{V,T} =
-\left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_T + \frac{RT}{V}
$$

```fortran
pressure: block
    real(pr) :: T, V, P, dPdV, dPdT, dPdn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 0.1_pr     ! Set volume to 0.1 L

    ! Calculate pressure and its derivatives
    call model%pressure(n, V, T, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)

    print *, "Pressure: ", P
    print *, "dPdV: ", dPdV
    print *, "dPdT: ", dPdT
    print *, "dPdn: ", dPdn

end block pressure
```

## Volume: \(V(n, P, T)\)
As you must know, given a temperature and pressure, the volume has not a unique
solution. In the case of a cubic EOS, there are three possible solutions for
the volume. For this reason, the volume is calculated iteratively. Provided
\(n\), \(P\), and \(T\) we fix \(n\) and
\(T\) and iterate over the volume until the pressure is the same as the
input pressure.

To learm how the initial guess for \(V\) is obtained, please refer to the book
"Thermodynamic Models: Fundamentals and Computational Aspects" by Michael L.
Michelsen and Jørgen M. Mollerup. Also, you will find how to calculate \(V\) derivatives in terms of \(P\) and \(T\) derivatives.

```fortran
volume: block
    real(pr) :: T, P, V

    T = 300.0_pr   ! Set temperature to 300 K
    P = 1.0_pr     ! Set pressure to 1 bar

    ! Calculate different volume roots
    call model%volume(n, P, T, V=V, root="liquid")
    print *, "Liquid volume: ", V

    call model%volume(n, P, T, V=V, root="vapor")
    print *, "Vapor volume: ", V

    call model%volume(n, P, T, V=V, root="stable")
    print *, "Stable volume root: ", V

end block volume
```

## Fugacity coefficients (V,T): \(\ln \phi_i (V,T)\)
Natural logarithm of fugacity coefficients specifing \(V\) and \(T\) are
calculated as follows:

$$
\ln \hat{\phi}_i = \frac{1}{RT} \left( \frac{\partial A^r}{\partial n_i} \right)_{V,T} - \ln Z
$$

Remember that the compressibility factor \(Z\) is calculated as:

$$
Z = \frac{PV}{nRT}
$$

Next, the derivatives:

$$
\left(\frac{\partial \ln \hat{\phi_i}}{\partial T} \right)_{P,n} = 
\frac{1}{RT} \left(\frac{\partial^2 A^r}{\partial T \partial n_i}\right)_V
- \frac{1}{RT^2} \left(\frac{\partial A^r}{\partial n_i}\right)_{V,T}
+ \frac{1}{T} + \frac{1}{RT} \frac{\left(\frac{\partial P}{\partial n_i} 
\right)_{V,T} \left(\frac{\partial P}{\partial T} \right)_{V,n}}{
\left(\frac{\partial P}{\partial V}\right)_{T,n}}
$$

$$
\left(\frac{\partial \ln \hat{\phi_i}}{\partial P} \right)_{T,n} = 
 \frac{1}{RT} \frac{\left(\frac{\partial P}{\partial n_i} 
\right)_{V,T}}{\left(\frac{\partial P}{\partial V}\right)_{T,n}} 
- \frac{1}{P}
$$

$$
\left(\frac{\partial \ln \hat{\phi_i}}{\partial n_j} \right)_{P,T} =
\frac{1}{n_T} + \frac{1}{RT} 
\left(\frac{\partial^2 A^r}{\partial n_i \partial n_j}
+ \frac{\left(\frac{\partial P}{\partial n_i} \right)_{V,T}
\left(\frac{\partial P}{\partial n_j} \right)_{V,T}}
{\left(\frac{\partial P}{\partial V} \right)_{T,n}}\right)
$$

With:

$$
n_T = \sum_i n_i
$$

```fortran
lnphi_vt: block
    real(pr) :: T, V, lnPhi(2), dlnPhidP(2), dlnPhidT(2), dlnPhidn(2,2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

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
    real(pr) :: T, V, lnPhi(2), P, dPdV, dPdT, dPdn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    call eos%lnphi_vt(&
      n, V, T, lnPhi=lnPhi, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
      )

end block lnphi_vt_p
```

## Fugacity coefficients (P,T): \(\ln \phi_i (P,T)\)

There is no need for a direct way of calculating natural logarithm the fugacity
coefficients specifing \(P\) and \(T\) since we can solve volume from those
specifications and then calculate the fugacity coefficients using the \(\ln
\phi_i (V,T)\) method. For that reason you can choose the root of the volume
that you want to use to calculate the fugacity coefficients in the same way as
the volume method.

```fortran

lnphi_pt: block
    real(pr) :: T, P, lnPhi(2), dlnPhidP(2), dlnPhidT(2), dlnPhidn(2,2)

    T = 300.0_pr   ! Set temperature to 300 K
    P = 1.0_pr     ! Set volume to 1 bar

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
    real(pr) :: T, P, V, lnPhi(2), dPdV, dPdT, dPdn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    P = 1.0_pr     ! Set volume to 1 bar

    call eos%lnphi_pt(&
      n, P, T, V=V, root_type="stable", lnPhi=lnPhi, &
      dPdV=dPdV, dPdT=dPdT, dPdn=dPdn &
      )

end block lnphi_pt_v

```

## Fugacity (V,T): \(\ln f_i (V,T)\)

Alternative way of calculating fugacity directly from the residual Helmholtz
free energy:

```fortran
fugacity_vt: block
    real(pr) :: T, V, lnf(2), dlnfdV(2), dlnfdT(2), dlnfdn(2,2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    call eos%lnfug_vt(&
      n, V, T, lnf, &
      dlnfdV, dlnfdT, dlnfdn, &
      )

end block fugacity_vt
```

## Residual entropy (V,T): \(S^r(n,V,T)\)

We start explaining the residual entropy because it is the easiest property to
calculate. Also, helps us to understand how to calculate the next properties.

The residual entropy is calculated as follows:

$$
S^r = - \left(\frac{\partial A^r}{\partial T} \right)_{V,n}
$$

And its derivatives:

$$
\left(\frac{\partial S^r}{\partial T} \right)_{V,n} =
- \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,n}
$$


$$
\left(\frac{\partial S^r}{\partial V} \right)_{T,n} =
- \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_n
$$

$$
\left(\frac{\partial S^r}{\partial n_i} \right)_{V,T} =
- \left(\frac{\partial^2 A^r}{\partial n_i \partial T} \right)_V
$$

```fortran
residual_entropy: block
    real(pr) :: T, V, S, dSdV, dSdT, dSdn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    ! Calculate residual entropy and its derivatives
    call eos%residual_entropy(n, V, T, Sr=Sr, SrT=SrT, SrV=SrV, Srn=Srn)

    print *, "Residual entropy: ", Sr
    print *, "SrV: ", SrV
    print *, "SrT: ", SrT
    print *, "Srn: ", Srn

end block residual_entropy
```

## Residual enthalpy (V,T): \(H^r(n,V,T)\)

The residual enthalpy is calculated as follows:

$$
H^r(n,V,T) = H^r(n,P,T) = A^r(n,V,T) + T S^r(n,V,T) + PV - nRT
$$

We know that pressure can be calculated as:

$$
P = - \left(\frac{\partial A^r}{\partial V} \right)_{T,n} + \frac{nRT}{V}
$$

Then, we can obtain:

$$
P - \frac{nRT}{V} = - \left(\frac{\partial A^r}{\partial V} \right)_{T,n}
$$

$$
PV - nRT = - V \left(\frac{\partial A^r}{\partial V} \right)_{T,n}
$$

And we also know how to calculate residual entropy as:

$$
S^r = - \left(\frac{\partial A^r}{\partial T} \right)_{V,n}
$$

Then, we can obtain the residual enthalpy as:

$$
H^r = A^r - T \left(\frac{\partial A^r}{\partial T} \right)_{V,n}
- V \left(\frac{\partial A^r}{\partial V} \right)_{T,n}
$$

And its derivatives:

$$
\left(\frac{\partial H^r}{\partial T} \right)_{V,n} =
- T \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,n}
- V \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_n
$$

$$
\left(\frac{\partial H^r}{\partial V} \right)_{T,n} =
- T \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_n
- V \left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,n}
$$

$$
\left(\frac{\partial H^r}{\partial n_i} \right)_{V,T} =
\left(\frac{\partial A^r}{\partial n_i} \right)_{V,T}
- T \left(\frac{\partial^2 A^r}{\partial T \partial n_i} \right)_V
- V \left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_T
$$

```fortran
residual_enthalpy: block
    real(pr) :: T, V, Hr, HrV, HrT, Hrn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    ! Calculate residual enthalpy and its derivatives
    call eos%residual_enthalpy(n, V, T, Hr=Hr, HrV=HrV, HrT=HrT, Hrn=Hrn)

    print *, "Residual enthalpy: ", Hr
    print *, "HrV: ", HrV
    print *, "HrT: ", HrT
    print *, "Hrn: ", Hrn

end block residual_enthalpy
```

## Residual Gibbs free energy (V,T): \(G^r(n,V,T)\)

The residual Gibbs free energy is calculated as follows:

$$
G^r(n,V,T) = A^r(n,V,T) + PV - nRT
$$

As with residual enthalpy we can easily deduce:

$$
G^r = A^r - V \left(\frac{\partial A^r}{\partial V} \right)_{T,n}
$$

And its derivatives:

$$
\left(\frac{\partial G^r}{\partial T} \right)_{V,n} =
\left(\frac{\partial A^r}{\partial T} \right)_{V,n}
- V \left(\frac{\partial^2 A^r}{\partial T \partial V} \right)_n
$$

$$
\left(\frac{\partial G^r}{\partial V} \right)_{T,n} =
- V \left(\frac{\partial^2 A^r}{\partial V^2} \right)_{T,n}
$$

$$
\left(\frac{\partial G^r}{\partial n_i} \right)_{V,T} =
\left(\frac{\partial A^r}{\partial n_i} \right)_{V,T}
- V \left(\frac{\partial^2 A^r}{\partial V \partial n_i} \right)_T
$$

```fortran
residual_gibbs: block
    real(pr) :: T, V, Gr, GrV, GrT, Grn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    ! Calculate residual Gibbs free energy and its derivatives
    call eos%residual_gibbs(n, V, T, Gr=Gr, GrV=GrV, GrT=GrT, Grn=Grn)

    print *, "Residual Gibbs free energy: ", Gr
    print *, "GrV: ", GrV
    print *, "GrT: ", GrT
    print *, "Grn: ", Grn

end block residual_gibbs
```

## Residual internal energy (V,T): \(U^r(n,V,T)\)

The residual internal energy is calculated as follows:

$$
U^r(n,V,T) = A^r(n,V,T) + T S^r(n,V,T)
$$

Therefore:

$$
U^r = A^r - T \left(\frac{\partial A^r}{\partial T} \right)_{V,n}
$$

And its derivatives:

$$
\left(\frac{\partial U^r}{\partial T} \right)_{V,n} =
- T \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,n}
$$

$$
\left(\frac{\partial U^r}{\partial V} \right)_{T,n} =
\left(\frac{\partial A^r}{\partial V} \right)_{T,n}
- T \left(\frac{\partial^2 A^r}{\partial V \partial T} \right)_n
$$

$$
\left(\frac{\partial U^r}{\partial n_i} \right)_{V,T} =
\left(\frac{\partial A^r}{\partial n_i} \right)_{V,T}
- T \left(\frac{\partial^2 A^r}{\partial T \partial n_i} \right)_V
$$

```fortran
residual_internal_energy: block
    real(pr) :: T, V, Ur, UrV, UrT, Urn(2)

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    ! Calculate residual internal energy and its derivatives
    call eos%internal_energyresidual(n, V, T, Ur=Ur, UrV=UrV, UrT=UrT, Urn=Urn)

    print *, "Residual internal energy: ", Ur
    print *, "UrV: ", UrV
    print *, "UrT: ", UrT
    print *, "Urn: ", Urn

end block residual_internal_energy
```

## Residual constant volume heat capacity (V,T): \(C_V^r(n,V,T)\)

The residual constant volume heat capacity is calculated as follows:

$$
C_V^r(n,V,T) = - T \left(\frac{\partial^2 A^r}{\partial T^2} \right)_{V,n}
$$

```fortran
residual_cv: block
    real(pr) :: T, V, Cv

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    ! Calculate residual constant volume heat capacity
    call eos%residual_cv(n, V, T, Cv=Cv)

    print *, "Residual constant volume heat capacity: ", Cv

end block residual_cv
```

## Residual constant pressure heat capacity (V,T): \(C_P^r(n,V,T)\)

The residual constant pressure heat capacity is calculated as follows:

$$
C_P^r(n,V,T) = C_V^r(n,V,T) - T \frac{\left(\frac{\partial P}
{\partial T} \right)^2_{V,n}}{\left(\frac{\partial P}
{\partial V} \right)_{T,n}} - nR
$$


```fortran
residual_cp: block
    real(pr) :: T, V, Cp

    T = 300.0_pr   ! Set temperature to 300 K
    V = 1.0_pr     ! Set volume to 1 liter

    ! Calculate residual constant pressure heat capacity
    call eos%residual_cp(n, V, T, Cp=Cp)

    print *, "Residual constant pressure heat capacity: ", Cp

end block residual_cp
```

# References
[1] 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
Fundamentals & computational aspects (2. ed). Tie-Line Publications.
