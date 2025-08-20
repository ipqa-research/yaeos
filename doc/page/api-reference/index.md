---
title: API Reference
---

[TOC]

# API Reference

This page provides a comprehensive overview of the yaeos Fortran API, organized by functionality.

## Core Modules

### Main Module

```fortran
use yaeos
```

The main module imports all essential components:
- [[yaeos__constants(module)]] - Constants and precision settings
- [[yaeos__models(module)]] - All thermodynamic models  
- [[yaeos__equilibria(module)]] - Phase equilibrium calculations

### Constants and Types

```fortran
use yaeos__constants

integer, parameter :: pr = real64        ! Working precision
real(pr), parameter :: R = 0.08314...    ! Gas constant [bar⋅L/(mol⋅K)]
type(KindEnum), parameter :: root_kinds  ! Phase type enumeration
```

## Models

### Base Model Types

#### ArModel
Abstract base class for residual Helmholtz energy models:

```fortran
class(ArModel), allocatable :: model

! Core methods
call model%residual_helmholtz(n, V, T, Ar, ArV, ArT, ...)  ! Fundamental method
call model%pressure(n, V, T, P, dPdV, dPdT, dPdn, ...)    ! Pressure calculation  
call model%lnphi_vt(n, V, T, lnphi, dlnphidv, ...)        ! Fugacity coefficients
call model%volume(n, P, T, V, root)                       ! Volume from P,T
```

#### GeModel  
Abstract base class for excess Gibbs energy models:

```fortran
class(GeModel), allocatable :: ge_model

call ge_model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, ...)  ! Fundamental method
call ge_model%activity_coefficient(n, T, lngamma, ...)     ! Activity coefficients
```

### Implemented Models

#### Cubic Equations of State

```fortran
! Peng-Robinson variants
model = PengRobinson76(tc, pc, w, [kij], [lij])           ! Original PR
model = PengRobinson78(tc, pc, w, [kij], [lij])           ! Modified alpha function

! Soave-Redlich-Kwong
model = SoaveRedlichKwong(tc, pc, w, [kij], [lij])

! RKPR (three-parameter)
model = RKPR(tc, pc, w, [kij], [lij])

! With custom alpha functions
model = PengRobinson76(tc, pc, w, alpha=my_alpha_func)

! With mixing rules
model = CubicEoS(tc, pc, w, alpha_func, mixing_rule)
```

#### Multifluid Models

```fortran
! GERG-2008 natural gas model
model = Gerg2008(components)
```

#### Excess Gibbs Models

```fortran
! NRTL
ge_model = NRTL(a_nrtl, b_nrtl, c_nrtl)

! UNIQUAC  
ge_model = UNIQUAC(q, r, a_uniquac, b_uniquac)

! UNIFAC (liquid-vapor)
ge_model = UNIFAC(groups, group_area, group_volume, group_interactions)
```

## Phase Equilibrium

### Flash Calculations

```fortran
use yaeos, only: flash, EquilibriumState

type(EquilibriumState) :: result

! PT Flash
call flash(model, n, T=T, P=P, equilibrium=result)

! TV Flash  
call flash(model, n, T=T, V=V, equilibrium=result)

! PV Flash
call flash(model, n, P=P, V=V, equilibrium=result)

! Access results
print *, "Number of phases:", result%phases
print *, "Phase fractions:", result%beta  
print *, "Compositions:", result%x
print *, "Volumes:", result%V
```

### Saturation Points

```fortran
use yaeos, only: saturation_pressure, saturation_temperature

real(pr) :: Psat, Tsat
type(EquilibriumState) :: sat_state

! Pure component saturation pressure
call saturation_pressure(model, T, Psat, equilibrium=sat_state)

! Pure component saturation temperature  
call saturation_temperature(model, P, Tsat, equilibrium=sat_state)

! Mixture bubble/dew points
call saturation_pressure(model, n, T, Psat, kind="bubble")
call saturation_pressure(model, n, T, Psat, kind="dew")
```

### Phase Envelopes

```fortran
use yaeos, only: pt_envelope_2ph, px_envelope_2ph

type(PTEnvelope2Ph) :: envelope
type(PXEnvelope2Ph) :: px_env

! PT envelope for binary mixture
call pt_envelope_2ph(model, n, envelope, kind="VLE")

! PX envelope at constant temperature
call px_envelope_2ph(model, n, T, px_env, kind="VLE")

! Access envelope points
print *, "Critical point:", envelope%critical_point
print *, "Pressures:", envelope%points%P
print *, "Temperatures:", envelope%points%T
```

### Critical Points

```fortran
use yaeos, only: critical_point, critical_line

real(pr) :: Tc, Pc, rhoc
type(CriticalLine) :: crit_line

! Pure component critical point
call critical_point(model, Tc, Pc, rhoc)

! Critical line for mixtures
call critical_line(model, z1, z2, crit_line)
```

## Utilities

### Mathematical Functions

```fortran
use yaeos__math

! Newton solver
call newton_solver(f, x0, x_solution, [jacobian], [tol])

! Linear system solver  
call solve_system(A, b, x)

! Continuation method
call continuation(f, x0, dx, points, [max_points])
```

### Consistency Testing

```fortran
use yaeos__consistency

! Check analytical derivatives
call check_ar_derivatives(model, n, V, T, [tolerance])
call check_ge_derivatives(ge_model, n, T, [tolerance])

! Numerical vs analytical comparison
call compare_derivatives(model, n, V, T)
```

### Automatic Differentiation

```fortran
use yaeos__adiff

! For developing new models
type(hyperdual) :: n_hd(nc), V_hd, T_hd, Ar_hd

! Convert to hyperdual numbers
n_hd = n
V_hd = V  
T_hd = T

! Calculate with automatic derivatives
call my_model_ar(n_hd, V_hd, T_hd, Ar_hd)
```

## Data Types

### Equilibrium Results

```fortran
type :: EquilibriumState
    integer :: phases                    ! Number of phases
    real(pr), allocatable :: beta(:)     ! Phase fractions
    real(pr), allocatable :: x(:,:)      ! Compositions [comp, phase]
    real(pr), allocatable :: V(:)        ! Phase volumes [L/mol]
    real(pr), allocatable :: rho(:)      ! Phase densities [mol/L]  
    real(pr) :: P                        ! Pressure [bar]
    real(pr) :: T                        ! Temperature [K]
end type
```

### Phase Envelope Points

```fortran
type :: PTEnvelopePoint2Ph
    real(pr) :: P, T                     ! Pressure [bar], Temperature [K]
    real(pr), allocatable :: x(:)        ! Liquid composition
    real(pr), allocatable :: y(:)        ! Vapor composition
    real(pr) :: beta                     ! Vapor fraction
end type

type :: PTEnvelope2Ph  
    type(PTEnvelopePoint2Ph), allocatable :: points(:)
    type(CriticalPoint), allocatable :: critical_points(:)
end type
```

### Substance Data

```fortran
type :: Substance
    character(len=:), allocatable :: name
    real(pr) :: tc      ! Critical temperature [K]
    real(pr) :: pc      ! Critical pressure [bar]  
    real(pr) :: w       ! Acentric factor
    real(pr) :: mw      ! Molecular weight [g/mol]
end type

type :: Substances
    type(Substance), allocatable :: substances(:)
end type
```

## Error Handling

yaeos uses optional arguments and return codes for error handling:

```fortran
! Most functions return error information through optional arguments
call flash(model, n, T=T, P=P, equilibrium=result, iters=iterations)

! Check convergence
if (iterations < 0) then
    print *, "Flash calculation failed to converge"
end if

! Some functions use ieee_arithmetic for NaN checking
use ieee_arithmetic, only: ieee_is_nan
if (ieee_is_nan(result)) then
    print *, "Calculation resulted in NaN"
end if
```

## Units

yaeos uses consistent units throughout:

- **Temperature**: Kelvin [K]
- **Pressure**: bar [bar]  
- **Volume**: Liters [L]
- **Molar amount**: moles [mol]
- **Energy**: bar⋅L (= 0.1 J = 100 mJ)
- **Gas constant**: R = 0.08314... bar⋅L/(mol⋅K)

## Memory Management

yaeos uses modern Fortran allocatable arrays:

```fortran
! Arrays are automatically allocated/deallocated
real(pr), allocatable :: n(:), tc(:), pc(:), w(:)

! Models use allocatable components
class(ArModel), allocatable :: model

! Clean up is automatic when variables go out of scope
```

## Performance Tips

1. **Reuse model objects**: Create once, use many times
2. **Provide good initial guesses**: Speeds up iterative calculations
3. **Use analytical derivatives**: Always faster than numerical
4. **Pre-allocate work arrays**: For repeated calculations
5. **Compiler optimization**: Use `-O3 -march=native` for production

## Example Usage Patterns

### Basic Property Calculation

```fortran
program basic_properties
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), V, T, P, lnphi(2)
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.5, 0.5]
    V = 2.0
    T = 350.0
    
    ! Calculate properties
    call model%pressure(n, V, T, P)
    call model%lnphi_vt(n, V, T, lnphi)
    
    print *, "P =", P, "bar"
    print *, "ln(φ) =", lnphi
end program
```

### Complete Phase Equilibrium Analysis

```fortran
program phase_analysis
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model  
    type(EquilibriumState) :: eq
    real(pr) :: n(2), T, P
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.3, 0.7]
    T = 350.0
    P = 15.0
    
    ! Flash calculation
    call flash(model, n, T=T, P=P, equilibrium=eq)
    
    ! Analyze results
    select case (eq%phases)
    case (1)
        print *, "Single phase"
    case (2) 
        print *, "Two phases:"
        print *, "  Vapor fraction:", eq%beta(2)
        print *, "  Liquid composition:", eq%x(:, 1)
        print *, "  Vapor composition:", eq%x(:, 2)
    end select
end program
```

## See Also

- [Installation Guide](installation.html)
- [Usage Examples](usage/index.html)  
- [Model Documentation](usage/eos/index.html)
- [Phase Equilibrium Guide](usage/phase_equilibrium/index.html)
- [Contributing Guide](contributing/index.html)