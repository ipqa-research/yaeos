---
title: Analytical derivatives
---

[TOC]

# Analytical Derivatives in yaeos

`yaeos` provides analytical derivatives for thermodynamic properties, which are essential for efficient phase equilibrium calculations and optimization tasks. Analytical derivatives are more accurate and computationally efficient than numerical differentiation.

## Why Analytical Derivatives?

1. **Accuracy**: No truncation errors from finite differences
2. **Performance**: Much faster than numerical differentiation
3. **Robustness**: Better convergence in iterative algorithms
4. **Consistency**: Exact relationships between derivatives

## Available Analytical Derivatives

### Residual Helmholtz Energy Derivatives

The fundamental derivatives that `yaeos` calculates are:

- \(\frac{\partial A^r}{\partial V}\) - First volume derivative
- \(\frac{\partial A^r}{\partial T}\) - First temperature derivative  
- \(\frac{\partial^2 A^r}{\partial V^2}\) - Second volume derivative
- \(\frac{\partial^2 A^r}{\partial T^2}\) - Second temperature derivative
- \(\frac{\partial^2 A^r}{\partial V \partial T}\) - Cross derivative
- \(\frac{\partial A^r}{\partial n_i}\) - Composition derivatives
- \(\frac{\partial^2 A^r}{\partial V \partial n_i}\) - Mixed derivatives
- \(\frac{\partial^2 A^r}{\partial T \partial n_i}\) - Mixed derivatives
- \(\frac{\partial^2 A^r}{\partial n_i \partial n_j}\) - Second composition derivatives

### Property Derivatives

From the Helmholtz derivatives, all other thermodynamic property derivatives are calculated:

```fortran
use yaeos

class(ArModel), allocatable :: model
real(pr) :: n(2), V, T, P
real(pr) :: dPdV, dPdT, dPdn(2)

! Calculate pressure with derivatives
call model%pressure(n, V, T, P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)

print *, "Pressure:", P, "bar"
print *, "dP/dV:", dPdV, "bar⋅mol/L"  
print *, "dP/dT:", dPdT, "bar/K"
print *, "dP/dn:", dPdn, "bar/mol"
```

## Using Derivatives in Phase Equilibrium

### Fugacity Coefficients

Fugacity coefficients and their derivatives are automatically calculated:

```fortran
real(pr) :: lnphi(2), dlnphidT(2), dlnphidP(2), dlnphidn(2,2)

! Get fugacity coefficients with derivatives
call model%lnphi_pt(n, P, T, lnphi, dlnPhidP=dlnphidP, &
                    dlnPhidT=dlnphidT, dlnPhidn=dlnphidn)
```

### Phase Equilibrium Jacobian

For faster flash calculations, analytical Jacobians are used:

```fortran
! Flash calculations use analytical derivatives internally
call flash(model, n, T=T, P=P, equilibrium=result)
! The algorithm uses analytical derivatives for:
! - Fugacity coefficients
! - Activity coefficients  
! - Pressure derivatives
```

## Implementing Analytical Derivatives in Custom Models

When creating custom models, you must provide analytical derivatives:

```fortran
subroutine residual_helmholtz(self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, &
                              Arn, ArVn, ArTn, Arn2)
    class(MyCustomModel), intent(in) :: self
    real(pr), intent(in) :: n(:), v, t
    real(pr), intent(out) :: Ar
    real(pr), optional, intent(out) :: ArV, ArT, ArTV, ArV2, ArT2
    real(pr), optional, intent(out) :: Arn(:), ArVn(:), ArTn(:)
    real(pr), optional, intent(out) :: Arn2(:,:)
    
    ! Calculate Ar and requested derivatives
    if (present(ArV)) ArV = ! dAr/dV
    if (present(ArT)) ArT = ! dAr/dT  
    if (present(ArV2)) ArV2 = ! d²Ar/dV²
    if (present(ArT2)) ArT2 = ! d²Ar/dT²
    if (present(ArTV)) ArTV = ! d²Ar/dVdT
    if (present(Arn)) Arn = ! dAr/dn
    if (present(ArVn)) ArVn = ! d²Ar/dVdn
    if (present(ArTn)) ArTn = ! d²Ar/dTdn  
    if (present(Arn2)) Arn2 = ! d²Ar/dn²
end subroutine
```

## Automatic Differentiation Support

`yaeos` also supports automatic differentiation for developing new models:

```fortran
use yaeos__adiff

! Use automatic differentiation to generate derivatives
! during model development, then implement analytical versions
```

## Derivative Consistency Checking

`yaeos` provides tools to verify your analytical derivatives:

```fortran
use yaeos__consistency

! Check if your model's derivatives are consistent
call check_derivatives_consistency(model, n, V, T)
```

## Examples

### Calculating All Pressure Derivatives

```fortran
program pressure_derivatives
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), V, T, P
    real(pr) :: dPdV, dPdT, dPdn(2)
    real(pr) :: d2PdV2, d2PdT2, d2PdVdT
    real(pr) :: tc(2), pc(2), w(2)
    
    ! Setup binary system
    tc = [305.32, 425.12]
    pc = [48.72, 37.96] 
    w = [0.0995, 0.2002]
    
    model = PengRobinson76(tc, pc, w)
    
    n = [0.5, 0.5]
    V = 2.0  ! L/mol
    T = 350.0  ! K
    
    ! Calculate all pressure derivatives
    call model%pressure(n, V, T, P, &
                       dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)
    
    print *, "Pressure:", P
    print *, "First derivatives:"
    print *, "  dP/dV =", dPdV  
    print *, "  dP/dT =", dPdT
    print *, "  dP/dn =", dPdn
end program
```

### Speed Comparison

```fortran
! Analytical derivatives are much faster:
! Analytical: ~1-2 μs per evaluation
! Numerical:  ~10-50 μs per evaluation (depending on method)

! This difference becomes critical in:
! - Flash calculations (hundreds of evaluations)  
! - Optimization routines (thousands of evaluations)
! - Phase envelope tracing (continuous calculations)
```

## Best Practices

1. **Always use analytical derivatives** when available
2. **Verify derivatives** during model development using consistency checks
3. **Profile your code** to ensure derivatives aren't the bottleneck
4. **Use automatic differentiation** during development, then implement analytical versions for production

## Mathematical Background

The analytical derivatives in `yaeos` are based on:

1. **Chain rule applications** for composite functions
2. **Implicit differentiation** for constraint equations  
3. **Thermodynamic relations** (Maxwell relations, etc.)
4. **Exact differential calculus** for consistent derivative networks

## See Also

- [[yaeos__consistency(module)]] - Derivative consistency checking
- [[yaeos__adiff(module)]] - Automatic differentiation support
- [[newmodels]] - Creating models with analytical derivatives
