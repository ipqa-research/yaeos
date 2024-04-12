---
title: Adding your own models
---

Most of thermodynamic properties calculated in `yaeos` heavily depend on
different kind of models and their respective derivatives.
Since obtaining the derivatives of complex models can be a tedious and
error-prone task. We provide two different ways of getting them automatically
(in some cases with some performance-cost), but there is also the possibility
of using anallyitical obtained expressions instead.

The calculation of thermodynamic properties heavily depends on 

On `yaeos` there are three different ways of adding your own model:W

[TOC]

# Residual Helmholtz models
Residual Helmholtz models are the basis to obtain the residual properties.

The main basis in `yaeos` to define a new object that extends the `abstract type` called `ArModel`. Which enforces the expected functionality of these
kind of models.

```fortran
use yaeos, only: ArModel

type, extends(ArModel) :: MyNewModel
end type
```

The definition of an `ArModel` expects that two procedures are defined:

- [[abs_residual_helmholtz(interface)]]: Procedure to calculate residual Helmholtz energy and it's derivatives
- [[abs_volume_initializer(interface)]]: Volume initializer to find a liquid root, given a pressure and temperature.

```fortran
use yaeos, only: ArModel

type, extends(ArModel) :: MyNewModel
contains
    procedure :: residual_helmholtz => an_Ar_implementation
    procedure :: volume_initializer => an_v0_implementation
end type
```

Satisfying those requirements, our model will be ready to make calculations!

```fortran
use yaeos, only: pressure
use my_model, only: MyNewModel

type(MyNewModel) :: model

! Assuming model parameters are set-up
call pressure(model, n, V, T, P)
```

## Using operator overloading automatic differentiation with hyperdual numbers

## Using source-to-source automatic differentiation with tapenade

## Using analytical obtained expressions
