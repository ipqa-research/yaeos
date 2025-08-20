---
title: Alpha functions
---

[TOC]

# Alpha Functions in Cubic Equations of State

Alpha functions (\(\alpha\)) are a crucial component of cubic equations of state. They describe the temperature dependence of the attractive parameter and allow cubic EoS to accurately represent the vapor pressure of pure components.

In `yaeos`, alpha functions are implemented as extensible derived types that can be easily substituted in cubic EoS models.

## Available Alpha Functions

### Soave Alpha Function

The Soave alpha function is one of the most widely used:

\[ \alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2 \]

where \(T_r = T/T_c\) is the reduced temperature and \(k\) is a parameter typically correlated with the acentric factor.

```fortran
use yaeos, only: AlphaSoave

type(AlphaSoave) :: alpha_func
real(pr) :: k(n) 

! For simple components, k can be calculated from acentric factor
k = 0.48 + 1.574*w - 0.176*w**2

alpha_func = AlphaSoave(k)
```

### RKPR Alpha Function

The RKPR (Redlich-Kwong-Peng-Robinson) alpha function:

\[ \alpha(T_r) = \left(\frac{3}{2 + T_r}\right)^k \]

```fortran
use yaeos, only: AlphaRKPR

type(AlphaRKPR) :: alpha_func
real(pr) :: k(n)

! k parameter needs to be fitted or correlated for each component
alpha_func = AlphaRKPR(k)
```

### Mathias-Copeman Alpha Function

A more flexible three-parameter alpha function:

```fortran
use yaeos, only: AlphaMathiasCopeman

type(AlphaMathiasCopeman) :: alpha_func
real(pr) :: c1(n), c2(n), c3(n)

! Parameters must be fitted to experimental data
alpha_func = AlphaMathiasCopeman(c1, c2, c3)
```

## Using Alpha Functions with Cubic EoS

Alpha functions are used when setting up cubic equations of state:

```fortran
use yaeos, only: CubicEoS, AlphaSoave, PengRobinson76

class(ArModel), allocatable :: model
type(AlphaSoave) :: alpha_func
real(pr) :: tc(2), pc(2), w(2), k(2)

! Component properties
tc = [190.6, 647.1]  ! Critical temperatures [K]
pc = [46.0, 220.6]   ! Critical pressures [bar]
w = [0.008, 0.344]   ! Acentric factors

! Create alpha function
k = 0.48 + 1.574*w - 0.176*w**2
alpha_func = AlphaSoave(k)

! Create EoS with custom alpha function
model = PengRobinson76(tc, pc, w, alpha=alpha_func)
```

## Creating Custom Alpha Functions

You can create custom alpha functions by extending the base `AlphaFunction` type:

```fortran
type, extends(AlphaFunction) :: MyCustomAlpha
    real(pr), allocatable :: param1(:)
    real(pr), allocatable :: param2(:)
contains
    procedure :: alpha => my_alpha_function
end type

subroutine my_alpha_function(self, Tr, a, dadt, dadt2)
    class(MyCustomAlpha), intent(in) :: self
    real(pr), intent(in) :: Tr(:)
    real(pr), intent(out) :: a(:), dadt(:), dadt2(:)
    
    ! Your custom alpha function implementation
    ! Must provide alpha, first and second temperature derivatives
end subroutine
```

## References

1. Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state. Chemical Engineering Science, 27(6), 1197-1203.
2. Mathias, P. M., & Copeman, T. W. (1983). Extension of the Peng-Robinson equation of state to complex mixtures. Fluid Phase Equilibria, 13, 91-108.