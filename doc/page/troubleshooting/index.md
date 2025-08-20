---
title: FAQ and Troubleshooting
---

[TOC]

# Frequently Asked Questions and Troubleshooting

This page addresses common questions and issues encountered when using yaeos.

## Installation and Setup

### Q: I'm getting "No Fortran compiler found" error

**A:** Ensure you have a modern Fortran compiler installed:

```bash
# Check if gfortran is available
gfortran --version

# Ubuntu/Debian
sudo apt install gfortran

# macOS with Homebrew  
brew install gfortran

# Windows with MSYS2
pacman -S mingw-w64-x86_64-gcc-fortran
```

### Q: fpm command not found

**A:** Install the Fortran Package Manager:

```bash
pip install --user fpm

# Or with conda
conda install -c conda-forge fpm

# Ensure Python user scripts are in PATH
export PATH=$PATH:~/.local/bin
```

### Q: Can't clone yaeos repository

**A:** Try these alternatives:

1. Use SSH instead of HTTPS: `git@github.com:ipqa-research/yaeos.git`
2. Download as ZIP file and use local path in `fpm.toml`
3. Use a specific release version instead of main branch

## Compilation Issues

### Q: Compilation fails with "Fortran 2008 features not supported"

**A:** Update your compiler to a more recent version:

- gfortran 8.0 or later
- ifort 18.0 or later  
- ifx (Intel's new compiler)

### Q: Out of memory during compilation

**A:** Try these solutions:

```bash
# Increase stack size
ulimit -s unlimited

# Use less aggressive optimization
export FFLAGS="-O1"

# Compile with more memory-efficient options
fpm build --flag "-fmax-stack-var-size=0"
```

### Q: Link errors with LAPACK/BLAS

**A:** Ensure linear algebra libraries are available:

```bash
# Ubuntu/Debian
sudo apt install liblapack-dev libblas-dev

# macOS with Homebrew
brew install openblas

# Or use Intel MKL
export LDFLAGS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
```

## Runtime Issues

### Q: "NaN" results in calculations

**Common causes and solutions:**

1. **Invalid input conditions**
   ```fortran
   ! Check temperature and pressure ranges
   if (T <= 0 .or. P <= 0) then
       print *, "Error: Invalid T or P"
   end if
   ```

2. **Poor model parameters**
   ```fortran
   ! Verify critical constants are reasonable
   if (any(tc <= 0) .or. any(pc <= 0)) then
       print *, "Error: Invalid critical properties"
   end if
   ```

3. **Convergence issues**
   ```fortran
   ! Provide better initial guesses
   call flash(model, n, T=T, P=P, equilibrium=eq, k0=k_initial)
   ```

### Q: Flash calculations not converging

**A:** Try these approaches:

1. **Check system conditions**
   ```fortran
   ! Ensure conditions are within reasonable range
   if (T > maxval(tc) .or. P > maxval(pc)) then
       print *, "Warning: Conditions near/above critical point"
   end if
   ```

2. **Use stability analysis**
   ```fortran
   use yaeos, only: min_tpd
   
   real(pr) :: tpd_result
   call min_tpd(model, n, P, T, tpd_result)
   
   if (tpd_result > 0) then
       print *, "System is stable (single phase)"
   end if
   ```

3. **Improve initial estimates**
   ```fortran
   ! Use Wilson K-factors as initial guess
   use yaeos, only: k_wilson
   
   real(pr) :: k_init(nc)
   k_init = k_wilson(tc, pc, w, T, P)
   
   call flash(model, n, T=T, P=P, equilibrium=eq, k0=k_init)
   ```

### Q: Pressure calculations giving unrealistic values

**A:** Common causes:

1. **Volume too small**: Molecules overlap
   ```fortran
   ! Check against excluded volume
   real(pr) :: b_total = sum(n * b_covolume)
   if (V < 1.1 * b_total) then
       print *, "Warning: Volume too small"
   end if
   ```

2. **Temperature too low**: Below practical range
   ```fortran
   if (T < 0.3 * minval(tc)) then
       print *, "Warning: Temperature very low"
   end if
   ```

3. **Model limitations**: EoS not suitable for conditions

## Model Setup Issues

### Q: How to choose the right equation of state?

**A:** Selection guidelines:

- **Hydrocarbons**: Peng-Robinson (PR76/PR78)
- **Polar compounds**: Consider NRTL + SRK mixing rules  
- **Natural gas**: GERG-2008
- **High pressure**: Peng-Robinson or RKPR
- **Near critical**: RKPR or GERG-2008

### Q: Binary interaction parameters (kij) missing

**A:** Options:

1. **Use default (kij = 0)**: For similar compounds
2. **Estimate from group contribution**: 
3. **Fit to experimental data**:
   ```fortran
   ! Example fitting structure
   real(pr) :: kij_matrix(nc, nc) = 0.0
   kij_matrix(1, 2) = 0.1  ! Fitted value
   kij_matrix(2, 1) = kij_matrix(1, 2)  ! Symmetric
   ```

### Q: Critical properties not available

**A:** Estimation methods:

1. **Group contribution methods**: Constantinou-Gani, Joback
2. **Corresponding states**: Use similar compounds
3. **Online databases**: NIST, DIPPR

## Performance Issues

### Q: Calculations are very slow

**A:** Optimization strategies:

1. **Compiler optimization**:
   ```bash
   fpm build --flag "-O3 -march=native -ffast-math"
   ```

2. **Avoid repeated model creation**:
   ```fortran
   ! Good: Create once, use many times
   model = PengRobinson76(tc, pc, w)
   do i = 1, n_points
       call model%pressure(n, V, T, P)
   end do
   
   ! Bad: Creating model in loop
   do i = 1, n_points
       model = PengRobinson76(tc, pc, w)  ! Don't do this!
       call model%pressure(n, V, T, P)
   end do
   ```

3. **Use analytical derivatives**: Always faster than numerical

4. **Profile your code**:
   ```bash
   # Compile with profiling
   fpm build --flag "-pg"
   ./your_program
   gprof your_program gmon.out > profile.txt
   ```

### Q: Memory usage is very high

**A:** Memory optimization:

1. **Use appropriate array sizes**: Don't over-allocate
2. **Deallocate when done**:
   ```fortran
   deallocate(large_array)
   ```
3. **Check for memory leaks**: Use valgrind or similar tools

## Numerical Issues

### Q: Results change with compiler/optimization

**A:** Possible causes:

1. **Floating-point precision**: Use consistent precision
   ```fortran
   use yaeos__constants, only: pr
   real(pr) :: all_variables  ! Use same precision everywhere
   ```

2. **Optimization issues**: Reduce optimization level
   ```bash
   fpm build --flag "-O1 -g"  # Lower optimization + debug info
   ```

3. **Mathematical instability**: Check algorithm convergence

### Q: Different results on different machines

**A:** Ensure consistent:

1. **Compiler versions**: Use same compiler
2. **Library versions**: Same BLAS/LAPACK versions  
3. **Input precision**: Use same number format
4. **Random seeds**: If using random numbers

## Advanced Troubleshooting

### Debug Mode

Enable detailed debugging:

```fortran
! In your code, add debug prints
print *, "Debug: T=", T, "P=", P, "n=", n

! Use ieee_arithmetic to check for special values
use ieee_arithmetic
if (ieee_is_nan(result)) print *, "NaN detected!"
if (.not. ieee_is_finite(result)) print *, "Infinite value!"
```

### Consistency Testing

Verify your models and derivatives:

```fortran
use yaeos__consistency

! Test analytical derivatives
call check_ar_derivatives(model, n, V, T)

! Compare analytical vs numerical
call compare_derivatives(model, n, V, T)
```

### Minimal Example for Bug Reports

When reporting issues, provide a minimal example:

```fortran
program minimal_example
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(2), V, T, P
    
    ! Simplest possible case that shows the problem
    model = PengRobinson76([305.32, 425.12], [48.72, 37.96], [0.0995, 0.2002])
    n = [0.5, 0.5]
    V = 2.0
    T = 350.0
    
    call model%pressure(n, V, T, P)
    print *, "Result:", P
    
    ! Include what you expected vs what you got
    print *, "Expected approximately: XXX"
end program
```

## Getting Help

If you still have issues:

1. **Search existing issues**: [GitHub Issues](https://github.com/ipqa-research/yaeos/issues)
2. **Create new issue** with:
   - Operating system and compiler version
   - Complete minimal example
   - Expected vs actual behavior
   - Error messages (if any)
3. **Include system info**:
   ```bash
   gfortran --version
   fpm --version
   uname -a
   ```

## Contributing Solutions

Found a solution to a common problem? Consider:

1. **Update this FAQ**: Submit a pull request
2. **Improve documentation**: Add examples to relevant sections
3. **Create examples**: Add to the examples directory
4. **Write tests**: Prevent regression of issues

## Performance Benchmarks

Typical performance expectations:

- **Flash calculation**: 1-10 ms per point
- **Pressure calculation**: 1-10 μs per point  
- **Derivative calculation**: 2-20 μs per point
- **Phase envelope**: 1-60 seconds per envelope

If your performance is significantly worse, there may be an optimization opportunity.