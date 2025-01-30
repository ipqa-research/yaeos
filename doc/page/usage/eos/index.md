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

The modern wat of defining an EOS is through a mathematical model that express
the residuel Helmholtz free energy as a function of the mole number, volume and
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
A^r(n, V, T) = -\int_{\infty}^{V} \left(P - \frac{nRT}{V} \right) dV
$$

This summary is based on the emblematic book Thermodynamic Models: Fundamentals
and Computational Aspects by Michael L. Michelsen and Jørgen M. Mollerup. This
book is a must-read for anyone interested in the subject. To honor the masters
of thermodynamics, we will use their notation, which is as follows:

$$
\frac{A^r(n, V, T)}{RT} = F
$$

\(F\) is the reduced residual Helmholtz free energy.

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

$$
P = -RT \left(\frac{\partial F}{\partial V} \right)_{T,n} 
+ \frac{nRT}{V}
$$

$$
\left(\frac{\partial P}{\partial V} \right)_{T,n} = 
-RT \left(\frac{\partial^2 F}{\partial V^2} \right)_{T,n} - \frac{nRT}{V^2}
$$

$$
\left(\frac{\partial P}{\partial T} \right)_{V,n} =
-RT \left(\frac{\partial^2 F}{\partial V \partial T} \right)_n + \frac{nR}{V}
$$

$$
\left(\frac{\partial P}{\partial n_i} \right)_{V,T} =
-RT \left(\frac{\partial^2 F}{\partial V \partial n_i} \right)_T + \frac{RT}{V}
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