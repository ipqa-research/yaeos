---
title: Using yaeos
ordered_subpage: eos
                 excessmodels
                 phase_equilibrium
                 newmodels
---

[TOC]

# Getting started

## Fortran
Maybe you've heard of Fortran as that old and cryptic language that everyone is
afraid of. Well, not anymore! Fortran is really easy to understand and has been
updated a lot in the recent decades. There is a fairly direct guide on the
[fortran-lang site](https://fortran-lang.org/learn/)

## Getting yaeos
`yaeos` is a Fortran library intended to be used as a
[`fpm`](fpm.fortran-lang.org) package (fpm: Fortran Package Manager), `fpm` can
be easily easily obtained with the Python package manager `pip` with a simple:

```bash
pip install --user fpm
```

With `fpm` installed you can create a new Fortran project by running:

```bash
fpm new your_project_name
```

A new directory with the name of your project will be created. 

You can include `yaeos` in your `fpm` project by adding it as a dependency on
your `fpm.toml` file by adding this:

```toml
[dependencies]
stdlib="*"
yaeos = {git="https://github.com/ipqa-research/yaeos"}
```

Or maybe you want a specific version:

```toml
[dependencies]
stdlib="*"
yaeos = {git="https://github.com/ipqa-research/yaeos", tag="2.0.0"}
```


# Setting up a model
On `yaeos` there is a series of models implemented, right now we include
Residual Helmholtz energy models (like Cubic Equations of State), but plan on
extening to a broader variety. Moreover, `yaeos` implements excess Gibbs energy
models like NRTL, UNIQUAC, liquid-vapor UNIFAC, and more.

In this example we'll show how a model in `yaeos` can be used. We'll take the
`Peng-Robinson` equation of state as an example, but all the implemented models
can be seen at [[yaeos__models(module)]]. Inside your `app/main.f90` file use

```fortran
program main
    use yaeos

    ! Set the variable `model` as a generic `ArModel`
    class(ArModel), allocatable :: model

    ! Set the variables that we're going to use as variable length arrays
    real(pr), allocatable :: n(:), tc(:), pc(:), w(:)

    n = [0.3, 0.7]    ! Number of moles of each component [mol]
    tc = [190, 310]   ! Critical temperatures [K]
    pc = [14, 30]     ! Critical pressures [bar]
    w = [0.001, 0.03] ! Acentric factors [-]

    ! Now we set our model as the PengRobinson76 equation of state.
    model = PengRobinson76(tc, pc, w)

end program
```

And then it's all set, now we've set the `model` variable to use in our
calculations. The parameter `pr` is the precision used by `yaeos`, with
`real64` as its value.

# Calculating thermodynamic properties
Some thermodynamic properties can be calculated with `yaeos` models, and we're
adding more! In this example we'll calculate a PV isotherm from our previously
defined model. For the sake of simplicity all the next code blocks are assumed
to be extensions of the previous one, before the `end program` sentence.

```fortran
pv_isotherm: block
    real(pr) :: v, t, p      ! Thermodynamic variables
    real(pr) :: v0, vf, dv   ! End and start volumes
    integer :: i, n_points   ! iteration variable and how many points to calc

    t = 300.0_pr   ! Set temperature to 300 K
    n_points = 10  ! Set the number of points to calculate

    ! Set initial and final volume points
    v0 = 0.1
    vf = 10
    dv = (vf - v0)/(n_points + 1)

    ! Build isotherm
    print *, "Volume [L], Pressure [bar]"

    do i = 0, n_points + 1
        ! Set new volume point
        v = v0 + i*dv

        ! Calculate pressure
        call model%pressure(n, v, t, p) ! <- Pressure is stored in `p`

        ! Print results as a CSV table
        print *, v, ",", p
    end do
end block pv_isotherm
```

You will obtain a table with the data:

| Volume [L] | Pressure [bar] |
|------------|---------------|
| 0.100      | 456.048       |
| 1.000      | 18.936        |
| 1.900      | 11.310        |
| 2.800      | 8.044         |
| 3.700      | 6.238         |
| 4.600      | 5.093         |
| 5.500      | 4.303         |
| 6.400      | 3.725         |
| 7.300      | 3.284         |
| 8.200      | 2.936         |
| 9.100      | 2.655         |
| 10.000     | 2.423         |


Also some useful derivatives are available when calculating each property, they
can be easily accessed as optional arguments of the routine. For example, to
obtain the derivative of pressure with respect to volume the line that
calculates pressure should be changed to:

```fortran        
call model%pressure(n, v, t, p, dPdV=dPdV) ! Calculate pressure and dPdV
```

In the next sections of this documentation are listed all the available
properties and derivatives that can be calculated with `yaeos` models. Also,
in further sections we'll show how to calculate phase equilibriums and how to
create new models.

# Units
If you were paying attention to the previous examples, you may have noticed
that the units of yaeos are defined according to the ideal gas constant R with
units:

\(R = 0.08314462618 \frac{bar \, L}{K \, mol}\)

Because of that, pressure must be specified as bar, volume as liters,
temperature as Kelvin and number of moles as moles. The returns of the
properties will be in the same units as the inputs. Energetic properties as
enthalpy, entropy, etc will have unit of \([bar \, L]\), which equals to
centijoules \([cJ]\).