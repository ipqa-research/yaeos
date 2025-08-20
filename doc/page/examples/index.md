---
title: Examples
---

[TOC]

# yaeos Examples and Tutorials

This page provides comprehensive examples demonstrating various yaeos capabilities.

## Getting Started Examples

### Basic Model Setup

```fortran
program basic_setup
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: tc(2), pc(2), w(2)
    
    ! Binary system: ethane + n-butane
    tc = [305.32, 425.12]  ! Critical temperatures [K]
    pc = [48.72, 37.96]    ! Critical pressures [bar]
    w = [0.0995, 0.2002]   ! Acentric factors
    
    ! Create Peng-Robinson model
    model = PengRobinson76(tc, pc, w)
    
    print *, "Model created:", model%name
end program
```

### Simple Property Calculation

```fortran
program property_calculation
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), V, T, P, lnphi(2)
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    
    ! Conditions
    n = [0.5, 0.5]  ! Composition [mol]
    V = 2.0         ! Volume [L]
    T = 350.0       ! Temperature [K]
    
    ! Calculate pressure
    call model%pressure(n, V, T, P)
    print *, "Pressure:", P, "bar"
    
    ! Calculate fugacity coefficients
    call model%lnphi_vt(n, V, T, lnphi)
    print *, "ln(phi):", lnphi
end program
```

## Phase Equilibrium Examples

### PT Flash Calculation

```fortran
program pt_flash_example
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    type(EquilibriumState) :: eq
    real(pr) :: n(2), T, P
    
    ! Setup binary system
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    
    ! Conditions
    n = [0.3, 0.7]  ! Overall composition
    T = 350.0       ! Temperature [K]
    P = 15.0        ! Pressure [bar]
    
    ! Perform flash calculation
    call flash(model, n, T=T, P=P, equilibrium=eq)
    
    ! Display results
    print *, "Number of phases:", eq%phases
    
    if (eq%phases == 2) then
        print *, "Two-phase system detected:"
        print *, "  Vapor fraction:", eq%beta(2)
        print *, "  Liquid composition:", eq%x(:, 1)
        print *, "  Vapor composition:", eq%x(:, 2)
        print *, "  Liquid volume:", eq%V(1), "L/mol"
        print *, "  Vapor volume:", eq%V(2), "L/mol"
    else
        print *, "Single phase system"
    end if
end program
```

### Saturation Pressure Calculation

```fortran
program saturation_example
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: T, Psat
    type(EquilibriumState) :: sat_eq
    
    ! Pure component: n-butane
    model = PengRobinson76([425.12], [37.96], [0.2002])
    
    ! Calculate saturation pressure at 350 K
    T = 350.0
    call saturation_pressure(model, [1.0_pr], T, Psat, equilibrium=sat_eq)
    
    print *, "Saturation pressure at", T, "K:", Psat, "bar"
    print *, "Liquid density:", 1.0/sat_eq%V(1), "mol/L"
    print *, "Vapor density:", 1.0/sat_eq%V(2), "mol/L"
end program
```

### Phase Envelope

```fortran
program phase_envelope_example
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    type(PTEnvel2) :: envelope
    real(pr) :: z(2)
    integer :: i
    
    ! Binary system
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    z = [0.3, 0.7]
    
    ! Calculate PT envelope
    call pt_envelope_2ph(model, z, envelope, kind="VLE")
    
    ! Print envelope points
    print *, "T [K], P [bar]"
    do i = 1, size(envelope%points)
        print *, envelope%points(i)%T, ",", envelope%points(i)%P
    end do
    
    ! Print critical point if found
    if (allocated(envelope%critical_points)) then
        print *, "Critical point:"
        print *, "  Tc =", envelope%critical_points(1)%T, "K"
        print *, "  Pc =", envelope%critical_points(1)%P, "bar"
    end if
end program
```

## Advanced Examples

### Custom Mixing Rules

```fortran
program custom_mixing_rules
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: tc(2), pc(2), w(2), kij(2,2), lij(2,2)
    
    ! Component properties
    tc = [305.32, 425.12]
    pc = [48.72, 37.96] 
    w = [0.0995, 0.2002]
    
    ! Binary interaction parameters
    kij = reshape([0.0, 0.05, 0.05, 0.0], [2, 2])
    lij = kij * 0.5  ! Volume interaction parameters
    
    ! Create model with interaction parameters
    model = PengRobinson76(tc, pc, w, kij, lij)
    
    print *, "Model with custom mixing rules created"
end program
```

### Multiple Models Comparison

```fortran
program model_comparison
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: pr76, pr78, srk
    real(pr) :: n(2), V, T, P_pr76, P_pr78, P_srk
    real(pr) :: tc(2), pc(2), w(2)
    
    ! Setup
    tc = [305.32, 425.12]
    pc = [48.72, 37.96]
    w = [0.0995, 0.2002]
    
    ! Create different models
    pr76 = PengRobinson76(tc, pc, w)
    pr78 = PengRobinson78(tc, pc, w)
    srk = SoaveRedlichKwong(tc, pc, w)
    
    ! Common conditions
    n = [0.5, 0.5]
    V = 2.0
    T = 350.0
    
    ! Compare pressures
    call pr76%pressure(n, V, T, P_pr76)
    call pr78%pressure(n, V, T, P_pr78)
    call srk%pressure(n, V, T, P_srk)
    
    print *, "Pressure comparison at V=", V, "L, T=", T, "K:"
    print *, "  PR76:", P_pr76, "bar"
    print *, "  PR78:", P_pr78, "bar"
    print *, "  SRK: ", P_srk, "bar"
end program
```

### Excess Gibbs Models

```fortran
program excess_gibbs_example
    use yaeos
    implicit none
    
    class(GeModel), allocatable :: nrtl_model
    real(pr) :: n(2), T, lngamma(2)
    real(pr) :: a_nrtl(2,2), b_nrtl(2,2), c_nrtl(2,2)
    
    ! NRTL parameters for water + ethanol
    ! (These are example parameters)
    a_nrtl = reshape([0.0, 3.458, -0.801, 0.0], [2,2])
    b_nrtl = reshape([0.0, -586.1, 246.18, 0.0], [2,2])
    c_nrtl = reshape([0.0, 0.3, 0.3, 0.0], [2,2])
    
    ! Create NRTL model
    nrtl_model = NRTL(a_nrtl, b_nrtl, c_nrtl)
    
    ! Calculate activity coefficients
    n = [0.5, 0.5]
    T = 298.15  ! K
    
    call nrtl_model%activity_coefficient(n, T, lngamma)
    
    print *, "Activity coefficients:"
    print *, "  Component 1:", exp(lngamma(1))
    print *, "  Component 2:", exp(lngamma(2))
end program
```

## Practical Applications

### Vapor Pressure Curve

```fortran
program vapor_pressure_curve
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: T, Psat, T_min, T_max, dT
    integer :: i, n_points
    
    ! Pure water
    model = PengRobinson76([647.1], [220.6], [0.344])
    
    ! Temperature range
    T_min = 300.0   ! K
    T_max = 600.0   ! K
    n_points = 31
    dT = (T_max - T_min) / (n_points - 1)
    
    print *, "Temperature [K], Pressure [bar]"
    do i = 1, n_points
        T = T_min + (i-1) * dT
        
        if (T < 647.1) then  ! Below critical temperature
            call saturation_pressure(model, [1.0_pr], T, Psat)
            print *, T, ",", Psat
        end if
    end do
end program
```

### Density Calculation

```fortran
program density_calculation
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), P, T, V_liquid, V_vapor, rho_l, rho_v
    type(EquilibriumState) :: eq
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.3, 0.7]
    P = 15.0
    T = 350.0
    
    ! Flash to get phase volumes
    call flash(model, n, T=T, P=P, equilibrium=eq)
    
    if (eq%phases == 2) then
        V_liquid = eq%V(1)  ! L/mol
        V_vapor = eq%V(2)   ! L/mol
        
        ! Calculate densities
        rho_l = 1.0 / V_liquid  ! mol/L
        rho_v = 1.0 / V_vapor   ! mol/L
        
        print *, "Liquid density:", rho_l, "mol/L"
        print *, "Vapor density:", rho_v, "mol/L"
        print *, "Density ratio:", rho_l / rho_v
    end if
end program
```

### PV Isotherm

```fortran
program pv_isotherm
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), T, V, P, V_min, V_max, dV
    integer :: i, n_points
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.5, 0.5]
    T = 350.0
    
    ! Volume range
    V_min = 0.1   ! L/mol
    V_max = 10.0  ! L/mol
    n_points = 50
    dV = (V_max - V_min) / (n_points - 1)
    
    print *, "Volume [L/mol], Pressure [bar]"
    do i = 1, n_points
        V = V_min + (i-1) * dV
        
        call model%pressure(n, V, T, P)
        print *, V, ",", P
    end do
end program
```

## Debugging and Validation Examples

### Derivative Checking

```fortran
program check_derivatives
    use yaeos
    use yaeos__consistency
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), V, T
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.5, 0.5]
    V = 2.0
    T = 350.0
    
    ! Check analytical derivatives against numerical
    call check_ar_derivatives(model, n, V, T)
    
    print *, "Derivative consistency check completed"
end program
```

### Stability Analysis

```fortran
program stability_analysis
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), P, T, tpd_result
    
    ! Setup
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.3, 0.7]
    P = 15.0
    T = 350.0
    
    ! Perform stability analysis
    call min_tpd(model, n, P, T, tpd_result)
    
    if (tpd_result > 0) then
        print *, "System is stable (single phase)"
        print *, "TPD value:", tpd_result
    else
        print *, "System is unstable (phase split expected)"
        print *, "TPD value:", tpd_result
    end if
end program
```

## Real-World Applications

### Natural Gas Composition Analysis

```fortran
program natural_gas_analysis
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(4), T, P  ! CH4, C2H6, C3H8, CO2
    real(pr) :: tc(4), pc(4), w(4)
    type(EquilibriumState) :: eq
    
    ! Component properties
    tc = [190.6, 305.3, 369.8, 304.1]  ! K
    pc = [46.0, 48.7, 42.5, 73.8]      ! bar
    w = [0.011, 0.100, 0.152, 0.225]   ! -
    
    ! Create model
    model = PengRobinson76(tc, pc, w)
    
    ! Typical natural gas composition (mol%)
    n = [0.85, 0.10, 0.03, 0.02]  ! Normalized
    
    ! Pipeline conditions
    T = 288.15  ! K (15°C)
    P = 70.0    ! bar
    
    ! Calculate phase behavior
    call flash(model, n, T=T, P=P, equilibrium=eq)
    
    if (eq%phases == 1) then
        print *, "Single gas phase at pipeline conditions"
        print *, "Density:", 1.0/eq%V(1), "mol/L"
    else
        print *, "Warning: Two phases detected!"
        print *, "This could indicate condensate formation"
    end if
end program
```

### Refinery Application

```fortran
program distillation_column_reboiler
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(3), T, P  ! Light, medium, heavy components
    real(pr) :: tc(3), pc(3), w(3)
    type(EquilibriumState) :: eq
    real(pr) :: K_values(3), alpha_relative(3)
    
    ! Hypothetical C5/C6/C7 mixture
    tc = [469.7, 507.6, 540.2]  ! K
    pc = [33.7, 30.3, 27.4]     ! bar
    w = [0.252, 0.301, 0.350]   ! -
    
    model = PengRobinson76(tc, pc, w)
    
    ! Reboiler composition
    n = [0.2, 0.5, 0.3]
    
    ! Reboiler conditions
    T = 400.0  ! K
    P = 2.0    ! bar
    
    call flash(model, n, T=T, P=P, equilibrium=eq)
    
    if (eq%phases == 2) then
        ! Calculate K-values
        K_values = eq%x(:, 2) / eq%x(:, 1)  ! y/x
        
        ! Relative volatilities (relative to heaviest)
        alpha_relative = K_values / K_values(3)
        
        print *, "Reboiler analysis:"
        print *, "  Vapor fraction:", eq%beta(2)
        print *, "  K-values:", K_values
        print *, "  Relative volatilities:", alpha_relative
    end if
end program
```

## Performance Examples

### Batch Processing

```fortran
program batch_flash_calculations
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), T_range(100), P_range(100)
    type(EquilibriumState) :: eq
    integer :: i, j
    real(pr) :: start_time, end_time
    
    ! Setup model once
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.5, 0.5]
    
    ! Generate conditions
    do i = 1, 100
        T_range(i) = 300.0 + i * 2.0  ! 300-500 K
        P_range(i) = 1.0 + i * 0.5    ! 1-50 bar
    end do
    
    ! Time the calculations
    call cpu_time(start_time)
    
    do i = 1, 100
        do j = 1, 100
            call flash(model, n, T=T_range(i), P=P_range(j), equilibrium=eq)
        end do
    end do
    
    call cpu_time(end_time)
    
    print *, "Completed 10,000 flash calculations"
    print *, "Total time:", end_time - start_time, "seconds"
    print *, "Average time per flash:", (end_time - start_time) / 10000.0 * 1000.0, "ms"
end program
```

## See Also

- [Installation Guide](../installation/index.html)
- [API Reference](../api-reference/index.html)
- [Troubleshooting](../troubleshooting/index.html)
- [GitHub Examples Directory](https://github.com/ipqa-research/yaeos/tree/main/example)