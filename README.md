[![Fortran](https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat)](https://fortran-lang.org)
[![fpm](https://img.shields.io/badge/fpm-Fortran_package_manager-734f96)](https://fpm.fortran-lang.org)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://ipqa-research.github.io/yaeos/)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://github.com/ipqa-research/yaeos/blob/main/LICENSE)
[![CI](https://github.com/fedebenelli/yaeos/actions/workflows/CI.yml/badge.svg)](https://github.com/ipqa-research/yaeos/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ipqa-research/yaeos/graph/badge.svg?token=IDJYKV8XK6)](https://codecov.io/gh/ipqa-research/yaeos)

<p align="center">
<img src="media/logo.png" width="50%"/>
</p>

There are multiple open source equations of state libraries, like:

- [Clapeyron](https://github.com/ClapeyronThermo/Clapeyron.jl) `julia`
- [FeOs](https://github.com/feos-org/feos) `rust` with `Python` bindings
- [teqp](https://github.com/usnistgov/teqp) `C++` with `Python` bindings
- [thermo](https://github.com/CalebBell/thermo) `python`
- [thermopack](https://github.com/thermotools/thermopack) `Fortran` with `Python` bindings
- [CoolProp](https://github.com/CoolProp/CoolProp) `C++` with `Python` bindings

Here we are presenting yet another (still in development) one, that tackles the
same problem just, in another way. Mostly exploiting the readability and
extensibility of Modern Fortran for scientists to have an easy way to implement
new thermodynamic models without dealing with lower-level languages but still
getting decent performance.
And also this framework provides the possibility of using analytically obtained
derivatives so both options are easily available.

> This is an experimental work in progress and we recommend the before
> mentioned libraries if you are intending to use some of this in real work.
> Big part of the code comes from a refactoring process of older codes so
> not all parts are easily readable, yet.

We focus mainly on that the addition of a new thermodynamic model as easily as
possible. Also providing our models too!

## Available models
- CubicEoS
   - SoaveRedlichKwong
   - PengRobinson76
   - PengRobinson78
- ExcessGibbs models
   - NRTL
   - UNIFAC VLE

## Available properties
- Bulk Properties
   - Volume(n, P, T)
   - Pressure(n, V, T)
- Residual Properties
   - H^R(n, V, T)
   - S^R(n, V, T)
   - G^R(n, V, T)
   - Cv^R(n, V, T)
   - Cp^R(n, V, T)
- Phase-Equilibria
   - FlashPT, FlashVT
   - Saturation points (bubble, dew and liquid-liquid)
   - Phase Envelope PT (isopleths)

# A little taste of `yaeos`
A lot of users get the bad picture of Fortran being old and archaic since most
of the codes they've seen are written in ancient `F77`.

```fortran
use yaeos, only: PengRobinson76, ArModel

integer, parameter :: n=2   ! Number of components
real(8) :: V, T, P, dPdN(n) ! variables to calculate
class(ArModel), allocatable :: model ! Model

real(pr) :: z(n), tc(n), pc(n), w(n), kij(n, n), lij(n, n)

z = [0.3, 0.7]
tc = [190., 310.]
pc = [14., 30.]
w = [0.001, 0.03]
kij = reshape([0., 0.1, 0.1, 0.], [n,n])
lij = kij / 2 

model =  PengRobinson76(tc, pc, w, kij, lij)

V = 1
T = 150

call model%pressure(z, V, T, P)
print *, P

! Obtain derivatives adding them as optional arguments! 
call model%pressure(model, z, V, T, P, dPdN=dPdN)
print *, dPdN
```

Examples of code with simple applications showing the capabilities of `yaeos`
can be found at [example/tutorials](example/tutorials). Each example can be run
with:

```bash
 fpm run --example <example name here>
```

Not providing any example will show all the possible examples that can be run.

# How to install/run it

##  Dependencies
`yaeos` needs to have both `lapack` and `nlopt` libraries on your system.

### Debian/Ubuntu-like
```bash
sudo apt install libnlopt-dev libblas-dev liblapack-dev
```

## Installing `yaeos`
`yaeos` is intended to use as a [`fpm`](fpm.fortran-lang.org) package. `fpm`
is the Fortran Package Manager, which automates the compilation and running
process of Fortran libraries/programs.

You can either:

- Generate a new project that uses `yaeos` as a dependency with:

> ```bash
> fpm new my_project
> ```
> 
> In the `fpm.toml` file add:
> 
> ```toml
> [dependencies]
> yaeos = {git="https://github.com/ipqa-research/yaeos"}
> ```

- Clone this repository and just modify the executables in the `app` directory

> ```bash
> git clone https://github.com/ipqa-research/yaeos
> cd yaeos
> fpm run
> ```

## Developing with vscode
If your intention is either to develop for `yaeos` or to explore in more detail
the library with debugging. We provide some predefined defuaults to work with
`vscode`. You can add them to the cloned repository by running:

```bash
git clone https://github.com/ipqa-research/vscode-fortran .vscode
```

From the project main directory 

## Available examples
In this repository we provide a series of examples of the different things that
can be calculated with `yaeos`. The source codes for the examples can be seen
at the [example/tutorials](example/tutorials) directory.

All the examples can be run with

```bash
fpm run --example <example_name_here>
```

# Including new models with Automatic Differentiation.

## Hyperdual Numbers autodiff
We are using the `hyperdual` module developed by 
[Philipp Rehner](https://github.com/prehner) 
and [Gernot Bauer](https://github.com/g-bauer)

> The automatic differentiation API isn't fully optimized yet so performance is
> much slower than it should be.

A complete implementation of the PR76 Equation of State can me found in
`example/adiff/adiff_pr76.f90`. Or in the documentation pages.

## Tapenade-based autodiff
It is also possible to differentiate with `tapenade`. Examples can be seen
in the documentation pages or in [The tools directory](tools/tapenade_diff/)

# Documentation
The latest API documentation for the `main` branch can be found
[here](https://ipqa-research.github.io/yaeos). This was generated from the source
code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford). We're
working in extending it more.
