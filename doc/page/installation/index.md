---
title: Installation and Setup
---

[TOC]

# Installation and Setup

This guide covers how to install and set up `yaeos` for your Fortran projects.

## Prerequisites

### Fortran Compiler

`yaeos` requires a modern Fortran compiler that supports Fortran 2008 standard or later:

- **gfortran** 8.0 or later (recommended)
- **ifort** (Intel Fortran) 18.0 or later  
- **ifx** (Intel's new Fortran compiler)
- **flang** (LLVM Fortran)

### Fortran Package Manager (fpm)

`yaeos` is distributed as an `fpm` package. Install `fpm` using Python:

```bash
pip install --user fpm
```

Alternatively, you can install `fpm` using conda:

```bash
conda install -c conda-forge fpm
```

## Installation Methods

### Method 1: Using fpm (Recommended)

This is the easiest way to use `yaeos` in your projects.

#### 1. Create a new fpm project

```bash
fpm new my_project
cd my_project
```

#### 2. Add yaeos as a dependency

Edit your `fpm.toml` file and add `yaeos` as a dependency:

```toml
[dependencies]
yaeos = { git = "https://github.com/ipqa-research/yaeos" }
```

For a specific version:

```toml
[dependencies]
yaeos = { git = "https://github.com/ipqa-research/yaeos", tag = "v4.3.0" }
```

#### 3. Use yaeos in your code

Create your main program in `app/main.f90`:

```fortran
program main
    use yaeos
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: tc(2), pc(2), w(2)
    
    ! Setup Peng-Robinson EoS for ethane + n-butane
    tc = [305.32, 425.12]  ! Critical temperatures [K]
    pc = [48.72, 37.96]    ! Critical pressures [bar]
    w = [0.0995, 0.2002]   ! Acentric factors
    
    model = PengRobinson76(tc, pc, w)
    
    print *, "yaeos model successfully created!"
    print *, "Model name:", model%name
end program main
```

#### 4. Build and run

```bash
fpm run
```

### Method 2: Git Submodule

For more control over the yaeos version:

```bash
# In your project directory
git submodule add https://github.com/ipqa-research/yaeos.git external/yaeos

# Update fpm.toml
echo '[dependencies]' >> fpm.toml
echo 'yaeos = { path = "external/yaeos" }' >> fpm.toml
```

### Method 3: Manual Download

1. Download the source code from [GitHub](https://github.com/ipqa-research/yaeos)
2. Extract to your desired location
3. Add the path to your `fpm.toml`:

```toml
[dependencies]
yaeos = { path = "/path/to/yaeos" }
```

## Configuration

### Database Path

Some yaeos features use parameter databases. Set the database path:

```fortran
use yaeos__constants, only: database_path

! Set the path to your database directory
database_path = "/path/to/yaeos/database"
```

### Precision Settings

`yaeos` uses double precision by default. The precision is controlled by:

```fortran
use yaeos__constants, only: pr

! pr is set to real64 (double precision)
real(pr) :: my_variable  ! This will be double precision
```

## Verification

### Basic Test

Create a simple test to verify your installation:

```fortran
program test_yaeos
    use yaeos
    use yaeos__constants, only: R, pr
    implicit none
    
    class(ArModel), allocatable :: model
    real(pr) :: n(1), V, T, P
    
    ! Create a simple pure component model
    model = PengRobinson76([647.1_pr], [220.6_pr], [0.344_pr]) ! Water
    
    ! Test pressure calculation
    n = [1.0_pr]  ! 1 mol
    V = 1.0_pr    ! 1 L
    T = 373.15_pr ! 100°C
    
    call model%pressure(n, V, T, P)
    
    print *, "Pressure calculated:", P, "bar"
    print *, "Expected around:", n(1)*R*T/V, "bar (ideal gas)"
    
    if (P > 0 .and. P < 1000) then
        print *, "✓ yaeos installation successful!"
    else
        print *, "✗ Something went wrong"
    end if
end program test_yaeos
```

Save this as `test/test_installation.f90` and run:

```bash
fpm test test_installation
```

### Running Examples

`yaeos` comes with examples in the `example/` directory:

```bash
# Run a basic example
fpm run --example basic_usage

# Run all examples
fpm run --example
```

## Troubleshooting

### Common Issues

#### Compilation Errors

**Problem**: Compiler not found
```
Error: No Fortran compiler found
```

**Solution**: Ensure a Fortran compiler is installed and in your PATH:
```bash
# Check if gfortran is available
gfortran --version

# On Ubuntu/Debian
sudo apt install gfortran

# On macOS with Homebrew
brew install gfortran

# On Windows with MSYS2
pacman -S mingw-w64-x86_64-gcc-fortran
```

**Problem**: Modern Fortran features not supported
```
Error: Fortran 2008 features not supported
```

**Solution**: Update your compiler to a more recent version.

#### FPM Issues

**Problem**: fpm command not found
```
bash: fmp: command not found
```

**Solution**: 
1. Check if fpm is installed: `pip list | grep fpm`
2. Ensure Python user scripts are in PATH
3. Reinstall: `pip install --user fpm --force-reinstall`

#### Network/Download Issues

**Problem**: Cannot clone yaeos repository
```
fatal: unable to access 'https://github.com/ipqa-research/yaeos/': ...
```

**Solution**: 
1. Check internet connection
2. Try using SSH instead: `git@github.com:ipqa-research/yaeos.git`
3. Download as ZIP file and use local path

### Performance Optimization

#### Compiler Flags

For production use, optimize compilation:

```toml
# In fpm.toml
[build]
auto-executables = true
auto-tests = true

[fortran]
implicit-typing = false

# Add compiler-specific optimizations
[compiler.gfortran]
flags = ["-O3", "-march=native", "-ffast-math"]

[compiler.ifort] 
flags = ["-O3", "-xHost", "-fast"]
```

#### Memory Usage

For large systems, you might need to increase stack size:

```bash
# Bash/zsh
ulimit -s unlimited

# Or set in your program
export OMP_STACKSIZE=1G
```

## Next Steps

After successful installation:

1. Read the [Usage Guide](usage/index.html)
2. Explore [Examples](https://github.com/ipqa-research/yaeos/tree/main/example)
3. Check out [Equation of State documentation](usage/eos/index.html)
4. Learn about [Phase Equilibrium calculations](usage/phase_equilibrium/index.html)

## Development Installation

If you plan to contribute to yaeos:

```bash
# Clone the development repository
git clone https://github.com/ipqa-research/yaeos.git
cd yaeos

# Run tests to ensure everything works
fpm test

# Generate documentation
ford ford.md

# Build examples
fpm run --example --all
```

See the [Contributing Guide](contributing/index.html) for more details.

## Support

If you encounter issues:

1. Check the [FAQ section](#troubleshooting) above
2. Search existing [GitHub Issues](https://github.com/ipqa-research/yaeos/issues)
3. Create a new issue with:
   - Your operating system
   - Compiler version (`gfortran --version`)
   - FPM version (`fpm --version`)
   - Minimal example that reproduces the problem