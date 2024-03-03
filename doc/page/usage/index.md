---
title: Using yaeos
---

[TOC]

# Getting started

## Fortran
Maybe you've heard of Fortran as that old and cryptic language that everyone is
afraid of. Well, not anymore! Fortran is really easy to understand and has been
updated a lot in the recent decades. There is a fairly direct guide on the
[fortran-lang site](https://fortran-lang.org/learn/)

## Getting yaeos
`yaeos`is a Fortran library intended to be used as a
[`fpm`](fpm.fortran-lang.org) package, `fpm` can be easily easily obtained with
the Python package manager `pip` with a simple:

```bash
pip install --user fpm
```

With `fpm` installed you can create a new Fortran project by running:

```bash
fpm new your_project_name
```

A new directory with the name of your project will be created. 

You include `yaeos` in your  `fpm` project by adding it as a dependency on your
`fpm.toml` file by adding this:

```toml
[dependencies]
stdlib="*"
yaeos = {git="https://github.com/ipqa-research/yaeos"}
```

Or maybe you want a specific version:

```toml
[dependencies]
stdlib="*"
yaeos = {git="https://github.com/ipqa-research/yaeos", tag="0.1.0b2"}
```


# Setting up a model
On `yaeos` there is a series of models implemented, right now we include
Residual Helmholtz energy models (like Cubic Equations of State), but plan on
extening to a broader variety.

In this example we'll show how a model in `yaeos` can be used. We'll take
the `Peng-Robinson` equation of state as an example, but all the implemented
models can be seen at [models](../../module/yaeos_models.html). Inside
your `app/main.f90` file use

```fortran
program main
    use yaeos

    ! Set the variable `model` as a generic `ArModel`
    class(ArModel), allocatable :: model

    ! Set the the variables that we're going to use
    ! as variable lenght arrays
    real(pr), allocatable :: n(:), tc(:), pc(:), w(:)

    n = [0.3, 0.7]    ! Number of moles of each component
    tc = [190, 310]   ! Critical temperatures
    pc = [14, 30]     ! Critical pressures
    w = [0.001, 0.03] ! Acentric factors

    ! Now we set our model as the PengRobinson76
    ! Equation of state.
    model = PengRobinson76(tc, pc, w)

end program
```

And then it's all set, now we've set the `model` variable to use in our
calculations

# Calculating thermodynamic properties
Some thermodynamic properties can be calculated with `yaeos` models, and we're
adding more! In this example we'll calculate a PV isotherm from our previously
defined model. For the sake of simplicity all the next code blocks are assumed
to be extensions of the previous one, before the `end program` sentence.

```fortran
pv_isotherm: block
    real(pr) :: v, t, p   ! Thermodynamic variables
    real(pr) :: v0, vf, dv ! End and start volumes
    integer :: i, n_points ! iteration variable and how many points to calc

    v0 = 0.001
    vf = 10
    dv = (vf - v0)/n_points

    do i=1,n_points
        v = v0 + i*dv  ! Set new volume point

        call pressure(model, n, v, t, p) ! Calculate pressure

        print *, v, p
    end do
end block pv_isotherm
```

Also some useful derivatives are available when calculating each property, they
can be easily accessed as optional arguments of the routine. For example, to 
obtain the derivative of pressure with respect to volume the line that
calculates pressure should be changed to:

```fortran        
call pressure(model, n, v, t, p, dpdv=dpdv) ! Calculate pressure and dPdV
```

The available thermodynamic properties to calculate can be seen at the
[Thermoprops](../../module/yaeos_thermoprops.html) module.

# Phase equilibria
