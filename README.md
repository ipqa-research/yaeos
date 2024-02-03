[![GitHub](https://img.shields.io/badge/GitHub-fedebenelli-blue.svg?style=social&logo=github)](https://github.com/fedebenelli) [![Fortran](https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat)](https://fortran-lang.org) [![fpm](https://img.shields.io/badge/fpm-Fortran_package_manager-734f96)](https://fpm.fortran-lang.org) [![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://fedebenelli.github.io/yaeos/) [![License](https://img.shields.io/github/license/fedebenelli/yaeos?color=green)](https://github.com/fedebenelli/yaeos/blob/main/LICENSE) [![CI](https://github.com/fedebenelli/yaeos/actions/workflows/CI.yml/badge.svg)](https://github.com/fedebenelli/yaeos/actions/workflows/CI.yml)

# YA_EoS
There are multiple open source equation of state libraries, like:

- [Clapeyron](https://github.com/ClapeyronThermo/Clapeyron.jl) `julia`
- [FeOs](https://github.com/feos-org/feos) `rust` with `Python` bindings
- [teqp](https://github.com/usnistgov/teqp) `C++` with `Python` bindings
- [thermo](https://github.com/CalebBell/thermo) `python`
- [thermopack](https://github.com/thermotools/thermopack) `Fortran` with `Python` bindings
- [CoolProp](https://github.com/CoolProp/CoolProp) `C++` with `Python` bindings

Here we are presenting yet another (still in development) one, that tackles the
same problem just, in another way. Mostly exploiting readability and
extensiblity of Modern Fortran for scientists to have an easy way to implement
new thermodynamic models without dealing with lower level languages but still
getting decent performance. 
And also this framework provides the possiblility of using analitically obtained
derivatives so both options are easy available.

> This is an experimental work in progress and we recommend the before
> mentioned libraries if you are intending to use some of this in real work.
> Big part of the code comes from a refactoring process of older codes so
> not all parts are easily readable, yet.

We focus mainly in that the addition of a new thermodynamic model is as easy as
possible. Also providing our own models too!

For now we only include residual helmholtz model (like Cubic or Saft Equations
of State). But we'll be adding other models like $G^E$ (UNIFAC for example).

## Including new models with Automatic Differentiation.
We are using the `hyperdual` module developed by [Philipp Rehner](https://github.com/prehner) and [Gernot Bauer](https://github.com/g-bauer)

> The automatic differentiation API isn't fully optimized yet so performance is
> much slower than it should be.

```fortran
type, extends(ArModelAdiff) :: YourNewModel
      type(Substances) :: composition
      real(8) :: parameters(:)
    contains
      procedure :: Ar => arfun
      procedure :: get_v0 => v0
end type
```

```fortran
subroutine arfun(self, n, v, t, Ar)
   class(YourNewModel), intent(in) :: self
   type(hyperdual), intent(in) :: n(:) ! Number of moles
   type(hyperdual), intent(in) :: v ! Volume [L]
   type(hyperdual), intent(in) :: t ! Temperature [K]
   type(hyperdual), intent(out) :: ar_value ! Residual Helmholtz Energy

   ! A very complicated residual helmholtz function of a mixture
   Ar = sum(n) * v * t
end subroutine

function v0(self, n, p, t)
   class(YourNewModel), intent(in) :: self
   real(pr), intent(in) :: n(:) ! Number of moles
   real(pr), intent(in) :: p ! Pressure [bar]
   real(pr), intent(in) :: t ! Temperature [K]
   real(pr) :: v0

   v0 = self%parameters(3)
end function
```

A complete implementation of the PR76 Equation of State can me found in
`example/adiff/adiff_pr76.f90`.

All the thermodynamic properties can be found in `yaeos_thermoprops` and called
like:

```fortran
use yaeos_thermoprops, only: fugacity_vt
use my_new_model, only: YourNewModel
...
type(YourNewModel) :: eos
eos%parameters = [1, 2, 3]
call fugacity_vt(eos, n, v, t, lnfug=lnfug, dlnphidn=dlnphidn)
```

## Documentation
The latest API documentation for the `main` branch can be found
[here](https://fedebenelli.github.io/yaeos). This was generated from the source
code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford). We're
working in extending it more.
