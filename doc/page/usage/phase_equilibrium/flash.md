---
title: Flash calculations
---

[TOC]

# Flash Calculations

Flash calculations are fundamental to phase equilibrium computations. They determine how many phases exist at given conditions and the composition of each phase.

## Types of Flash Calculations

### PT Flash (Pressure-Temperature Flash)

The most common flash calculation where pressure, temperature, and overall composition are specified:

```fortran
use yaeos, only: flash, EquilibriumState

class(ArModel), allocatable :: model
type(EquilibriumState) :: result
real(pr) :: n(nc), T, P

! Setup your model and initial composition
n = [0.3, 0.7]  ! Overall composition [mol]
T = 350.0       ! Temperature [K]
P = 20.0        ! Pressure [bar]

! Perform PT flash
call flash(model, n, T=T, P=P, equilibrium=result)

! Check results
if (result%phases == 1) then
    print *, "Single phase"
else if (result%phases == 2) then
    print *, "Two phases detected"
    print *, "Vapor fraction:", result%beta(2)
    print *, "Liquid composition:", result%x(:, 1)
    print *, "Vapor composition:", result%x(:, 2)
end if
```

### TV Flash (Temperature-Volume Flash)

Flash calculation at specified temperature and volume:

```fortran
real(pr) :: V = 2.5  ! Total volume [L]

call flash(model, n, T=T, V=V, equilibrium=result)
```

### PV Flash (Pressure-Volume Flash)

Flash calculation at specified pressure and volume:

```fortran
call flash(model, n, P=P, V=V, equilibrium=result)
```

## Understanding Flash Results

The `EquilibriumState` type contains all the information about the equilibrium:

```fortran
type(EquilibriumState) :: eq

! Number of phases present
print *, "Number of phases:", eq%phases

! Phase fractions (beta values)
print *, "Phase fractions:", eq%beta

! Compositions of each phase
do i = 1, eq%phases
    print *, "Phase", i, "composition:", eq%x(:, i)
end do

! Densities and volumes
print *, "Phase densities [mol/L]:", eq%rho
print *, "Phase volumes [L]:", eq%V
```

## Multiphase Flash

For systems with more than two phases, use the multiphase flash:

```fortran
use yaeos, only: pt_mp_flash, MPEquilibriumState

type(MPEquilibriumState) :: mp_result

call pt_mp_flash(model, n, T, P, mp_result)

print *, "Number of phases found:", size(mp_result%compositions, 2)
```

## Flash Calculation Tips

### Initial Estimates

Good initial estimates improve convergence:

```fortran
! Provide initial guess for phase split
real(pr) :: beta_init(2) = [0.3, 0.7]  ! 30% liquid, 70% vapor

call flash(model, n, T=T, P=P, equilibrium=result, beta0=beta_init)
```

### Convergence Control

You can control the flash calculation tolerance:

```fortran
! More strict convergence
call flash(model, n, T=T, P=P, equilibrium=result, tol=1e-10)
```

### Handling Convergence Issues

If flash calculations fail to converge:

1. **Check your model setup**: Ensure critical properties and interaction parameters are reasonable
2. **Verify conditions**: Are you within the reasonable range for your system?
3. **Try different initial guesses**: Sometimes a better starting point helps
4. **Use stability analysis**: Check if the phase split is thermodynamically stable

```fortran
use yaeos, only: min_tpd

real(pr) :: tpd_result
logical :: is_stable

! Check stability of single phase
call min_tpd(model, n, P, T, tpd_result)
is_stable = tpd_result > 0
```

## Examples

### Binary System Flash

```fortran
program binary_flash
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    type(EquilibriumState) :: eq
    real(pr) :: n(2), T, P, tc(2), pc(2), w(2)
    
    ! Setup PR EoS for ethane + n-butane
    tc = [305.32, 425.12]    ! K
    pc = [48.72, 37.96]      ! bar  
    w = [0.0995, 0.2002]
    
    model = PengRobinson76(tc, pc, w)
    
    ! Flash calculation
    n = [0.5, 0.5]
    T = 350.0  ! K
    P = 15.0   ! bar
    
    call flash(model, n, T=T, P=P, equilibrium=eq)
    
    if (eq%phases == 2) then
        print *, "Vapor fraction:", eq%beta(2)
        print *, "Liquid mole fractions:", eq%x(:, 1) 
        print *, "Vapor mole fractions:", eq%x(:, 2)
    end if
end program
```

## See Also

- [[yaeos__equilibria_flash(module)]] - API documentation
- [[saturation_points]] - For pure component saturation calculations
- [[envelopes]] - For phase envelope calculations