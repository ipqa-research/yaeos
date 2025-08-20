---
title: Saturation Points
---

[TOC]

# Saturation Points Calculations

Saturation points are conditions where a pure component or mixture exists in equilibrium between two phases (typically liquid and vapor). These calculations are fundamental for understanding phase behavior.

## Pure Component Saturation

### Saturation Pressure

Calculate the saturation pressure at a given temperature:

```fortran
use yaeos, only: saturation_pressure

class(ArModel), allocatable :: model
real(pr) :: T, Psat
real(pr) :: tc, pc, w

! Setup pure component (e.g., n-butane)
tc = 425.12  ! Critical temperature [K]
pc = 37.96   ! Critical pressure [bar]
w = 0.2002   ! Acentric factor

model = PengRobinson76([tc], [pc], [w])

! Calculate saturation pressure at 350 K
T = 350.0
call saturation_pressure(model, T, Psat)

print *, "Saturation pressure at", T, "K:", Psat, "bar"
```

### Saturation Temperature

Calculate the saturation temperature at a given pressure:

```fortran
use yaeos, only: saturation_temperature

real(pr) :: P, Tsat

! Calculate saturation temperature at 10 bar
P = 10.0
call saturation_temperature(model, P, Tsat)

print *, "Saturation temperature at", P, "bar:", Tsat, "K"
```

## Mixture Saturation Points

For mixtures, saturation calculations are more complex as they involve phase equilibrium between phases of different compositions.

### Dew Point Calculation

The dew point is the temperature (or pressure) at which the first drop of liquid appears:

```fortran
use yaeos, only: saturation_temperature

real(pr) :: z(2), P, T_dew
type(EquilibriumState) :: eq

! Overall composition
z = [0.3, 0.7]
P = 15.0  ! bar

! Find dew point temperature
call saturation_temperature(model, P, T_dew, z, kind="dew")

print *, "Dew point temperature:", T_dew, "K"
```

### Bubble Point Calculation

The bubble point is the temperature (or pressure) at which the first bubble of vapor appears:

```fortran
real(pr) :: T_bubble

! Find bubble point temperature  
call saturation_temperature(model, P, T_bubble, z, kind="bubble")

print *, "Bubble point temperature:", T_bubble, "K"
```

## Advanced Saturation Calculations

### Pure Component Vapor Pressure Curve

Generate a complete vapor pressure curve:

```fortran
program vapor_pressure_curve
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: tc, pc, w, T, Psat
    integer :: i, n_points
    
    ! Setup water
    tc = 647.1   ! K
    pc = 220.6   ! bar
    w = 0.344
    
    model = PengRobinson76([tc], [pc], [w])
    
    n_points = 20
    
    print *, "Temperature [K], Pressure [bar]"
    do i = 1, n_points
        T = 300.0 + real(i-1) * (tc - 300.0) / real(n_points-1)
        
        if (T < tc) then
            call saturation_pressure(model, T, Psat)
            print *, T, ",", Psat
        end if
    end do
end program
```

### Saturation Properties

You can also calculate other properties along the saturation curve:

```fortran
use yaeos, only: saturation_pressure

real(pr) :: T, Psat, Vl, Vv, rhol, rhov
type(EquilibriumState) :: sat_state

T = 350.0
call saturation_pressure(model, T, Psat, equilibrium=sat_state)

! Get liquid and vapor volumes
Vl = sat_state%V(1)  ! Liquid volume [L/mol]
Vv = sat_state%V(2)  ! Vapor volume [L/mol]

! Calculate densities
rhol = 1.0/Vl  ! Liquid density [mol/L]
rhov = 1.0/Vv  ! Vapor density [mol/L]

print *, "At T =", T, "K:"
print *, "Saturation pressure:", Psat, "bar"
print *, "Liquid density:", rhol, "mol/L"  
print *, "Vapor density:", rhov, "mol/L"
```

## Critical Points

The critical point is where the saturation curve terminates:

```fortran
use yaeos, only: critical_point

real(pr) :: tc_calc, pc_calc, rhoc_calc

! Calculate critical point
call critical_point(model, tc_calc, pc_calc, rhoc_calc)

print *, "Critical temperature:", tc_calc, "K"
print *, "Critical pressure:", pc_calc, "bar"
print *, "Critical density:", rhoc_calc, "mol/L"
```

## Practical Tips

### Convergence Issues

If saturation calculations fail to converge:

1. **Check your model**: Ensure critical properties are realistic
2. **Stay within bounds**: Don't calculate saturation properties too close to the critical point
3. **Use good initial guesses**: Provide reasonable starting values

### Initialization

For pure components, Antoine equation parameters can provide good initial estimates:

```fortran
! Antoine equation: log10(Psat) = A - B/(T + C)
real(pr) :: A, B, C, P_antoine

! Example for n-butane
A = 4.35576
B = 1175.581
C = -2.071

P_antoine = 10**(A - B/(T + C))  ! mmHg
P_antoine = P_antoine * 0.00133322  ! Convert to bar
```

## Examples

### Complete Saturation Analysis

```fortran
program saturation_analysis
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: tc, pc, w
    real(pr) :: T, Psat, T_target, P_target
    type(EquilibriumState) :: eq
    
    ! Setup component (propane)
    tc = 369.89
    pc = 42.48
    w = 0.1521
    
    model = PengRobinson76([tc], [pc], [w])
    
    ! Calculate Psat at different temperatures
    print *, "Vapor Pressure Curve:"
    print *, "T [K], Psat [bar]"
    
    do T = 200.0, 360.0, 20.0
        call saturation_pressure(model, T, Psat)
        print *, T, ",", Psat
    end do
    
    ! Calculate Tsat at different pressures  
    print *, ""
    print *, "Saturation Temperature:"
    print *, "P [bar], Tsat [K]"
    
    do P_target = 1.0, 35.0, 5.0
        call saturation_temperature(model, P_target, T_target)
        print *, P_target, ",", T_target
    end do
end program
```

## See Also

- [[yaeos__equilibria_saturation_points(module)]] - API documentation
- [[flash]] - For general phase equilibrium calculations
- [[envelopes]] - For mixture phase boundaries