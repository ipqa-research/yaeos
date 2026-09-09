---
title: User documentation
ordered_subpage: usage
                 contributing
---

# Welcome to `yaeos` user documentation

<p align="center">
  <img src="../media/logo.svg" width="500">
</p>

`yaeos` is a Fortran library that handles thermodynamic Equations of
State-based calculations, mostly phase-equilibria related ones and properties
estimation ones.

This is a work-in-progress library (and documentation) so don't hesisate to
report any problem/requirement as an issue in our [GitHub
page](https://github.com/ipqa-research/yaeos/issues).

## Basic usage
The base object that represents most of `yaeos` functionality is the `ArModel`
object, which holds the basic interface of how a \( A_r \) model behaves. 

Since all the properties that `yaeos` calculates are based on residual
Helmholtz energy, once the object is set-up all the library functionality is
available.

```fortran
use yaeos

class(ArModel), allocatable :: model ! Variable that holds the model

! A setup function that returns a setted model
model = setup_model(<the properties that define a model>)

! Once the model is set up, the user has access to the properties
call model%pressure(n, V, T, P, dPdn=dPdn)
```