[![GitHub](https://img.shields.io/badge/GitHub-fedebenelli-blue.svg?style=social&logo=github)](https://github.com/fedebenelli)

[![Fortran](https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat)](https://fortran-lang.org)

[![fpm](https://img.shields.io/badge/fpm-Fortran_package_manager-734f96)](https://fpm.fortran-lang.org)

[![Tag](https://img.shields.io/github/v/tag/fedebenelli/fortime?color=blue&logo=github&style=flat)](https://github.com/fedebenelli/yaeos/releases)

[![Release](https://img.shields.io/github/release/fedebenelli/yaeos.svg)](https://github.com/fedebenelli/yaeos/releases/latest)

[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://fedebenelli.github.io/yaeos/)

[![License](https://img.shields.io/github/license/fedebenelli/yaeos?color=green)](https://github.com/fedebenelli/yaeos/blob/main/LICENSE)

[![CI](https://github.com/fedebenelli/yaeos/actions/workflows/CI.yml/badge.svg)](https://github.com/fedebenelli/yaeos/actions/workflows/CI.yml)


# YA_EoS
There are multiple equation of state libraries, like:

- [Clapeyron](https://github.com/ClapeyronThermo/Clapeyron.jl) `julia`
- [FeOs](https://github.com/feos-org/feos) `rust`
- [OSSEoS](https://hpp.uva.es/open-source-software-eos/) `MATLAB`
- [teqp](https://github.com/usnistgov/teqp) `C++`
- [thermo](https://github.com/CalebBell/thermo) `python`
- [thermopack](https://github.com/thermotools/thermopack) `Fortran`

Here we are presenting yet another (still in development) one, that tackles the
same problem just, in another way. Mostly exploiting readability and
extensiblity of Modern Fortran for scientists to have an easy way to implement
new thermodynamic models without dealing with lower level languages but still
getting decent performance. And also having the possiblility of using
analitically obtained derivatives so both options are available with just a 
switch.

> This is an experimental work in progress and we recommend the before
> mentioned libraries if you are intending to use some of this in real work.

We focus mainly in that the addition of a new thermodynamic model is as easy as
possible.
Also providing our own models too!

For now we only include residual helmholtz model (like Cubic or Saft Equations
of State). But we'll be adding other models like $G^E$ (UNIFAC for example).

```fortran
subroutine Ar(z, v, t, ar)
   use yaeos_adiff
   type(hyperdual), intent(in) :: z(:)
   type(hyperdual), intent(in) :: v
   type(hyperdual), intent(in) :: t
   type(hyperdual), intent(out) :: ar

   ! A very complicated residual helmholtz function of a mixture
   ar = sum(z) * v * t
end subroutine
```

After defining your own model you must just setup the Ar function and you're
done, the available routines are using your model:

```fortran
use yaeos_ar_models, only: set_ar_function
use yaeos_thermoprops, only: pressure
...

call set_ar_function(Ar) ! Set the before defined subroutine
call pressure(z, v, t, p, dp, dp2) ! Get the pressure and it's derivatives 
                                   ! at some thermodynamic state
```

## Documentation
The latest API documentation for the `main` branch can be found
[here](https://fedebenelli.github.io/yaeos). This was generated from the source
code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).
