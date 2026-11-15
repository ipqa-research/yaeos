---
title: Architecture and repository map
---

[TOC]

# Find your way around the codebase

The central design idea in `yaeos` is to keep a thermodynamic model separate
from the calculations that use it. An equilibrium algorithm needs pressure,
fugacity, or derivatives; it should not have to know whether they came from
Peng–Robinson, GERG2008, or another compatible energy model.

That separation is implemented through Fortran modules, abstract derived types,
and type-bound procedures. This page explains those layers and gives a
concrete route from a user call to the mathematical kernel.

## The repository in one table

Paths in this handbook are relative to the repository root.

| Path | Responsibility | Start here when… |
|------|----------------|-----------------|
| `src/yaeos.f90` | Top-level public imports | A name is missing from `use yaeos` |
| `src/constants.f90` | Precision, gas constant, selectors, shared settings | Units, kinds, or root conventions matter |
| `src/models/` | Model base types and implementations | Adding or diagnosing property calculations |
| `src/equilibria/` | Stability, saturation, flash, critical and boundary calculations | A phase calculation is wrong or fails |
| `src/math/` | Linear systems, nonlinear solvers, continuation, numerical helpers | Diagnosing an algorithm's numerical machinery |
| `src/adiff/` | Hyperdual arithmetic and derivative adapters | Implementing automatic differentiation |
| `src/consistency/` | Derivative and thermodynamic consistency checks | Validating a new energy function |
| `src/extra/` | Example/test fluid fixtures | Reusing a documented mixture |
| `src/legacy/` | Older global-state implementations | Maintaining an existing legacy caller |
| `c_interface/yaeos_c.f90` | Flat wrapper procedures and model-handle storage | Changing a language-facing call |
| `python/yaeos/` | Python classes, conversion, fitting and tools | Changing Python user behavior |
| `test/` and `python/tests/` | Native and Python regression tests | Proving behavior across layers |
| `example/` | Demonstration programs | Learning a call sequence or adding a runnable example |
| `tools/` | Development/generation utilities | Regenerating differentiated code |
| `doc/page/` and `python/docs/source/` | FORD guides and Sphinx documentation | Explaining a feature to users |

The actual compiled source sets differ by build. `fpm.toml` describes the
native package; `python/meson.build` adds the interface and extension build.
Finding a file in the repository does not establish that every target builds
it. See [workflow](workflow.html) for the separate build paths.

## Public modules are entry points, not separate model implementations

`src/yaeos.f90` imports broad API groups. `src/models/models.f90` gathers model
types and constructors; `src/equilibria/equilibria.f90` gathers equilibrium
routines and result types. These modules are often called **facades**: they
make common names convenient to import without duplicating the underlying
calculation.

For example, these are alternative import styles, not different flash
algorithms:

```fortran
use yaeos
```

```fortran
use yaeos__equilibria, only: flash
```

The directory hierarchy is reflected in names such as
`yaeos__models_ar_genericcubic`. These are ordinary modules. Do not confuse an
organizational “submodule” in prose with the Fortran language's specific
`submodule` construct.

When adding a public feature, check both its defining module's visibility and
the relevant facade's `use ..., only:` list. Implementing a procedure does not
automatically make it available through every facade or through Python.

## The model family tree

```text
Substances: component metadata
    stored in BaseModel%components

BaseModel
    +-- ArModel: residual Helmholtz energy A^r(n,V,T)
    |     +-- CubicEoS
    |     +-- ArModelAdiff -> hyperdual models, including Gerg2008
    |     +-- ArModelTapenade -> source-differentiated models, including PCSAFT
    |
    +-- GeModel: excess Gibbs energy G^E(n,T)
          +-- analytical activity models, such as NRTL and UNIQUAC
          +-- hyperdual/source-differentiation adapter families
```

The tree shows implementation families, not a list of every concrete type.
`Substances` is a container of arrays in component order, not a parent class
of `BaseModel`.

### Why use an energy-function contract?

An `ArModel` supplies residual Helmholtz energy and requested derivatives.
Shared methods use them to calculate pressure, fugacity coefficients, residual
enthalpy, and other properties. A new compatible model can therefore reuse
property relations and equilibrium algorithms instead of independently
implementing each property.

`GeModel` follows the same idea for excess Gibbs energy. Its input does not
include volume or pressure in the base energy interface. It is not an EOS
implementation with a few arguments omitted: it describes a different
thermodynamic reference and different capabilities.

The detailed [Fortran contract guide](fortran.html) explains which methods are
deferred, how optional derivatives work, and why extensive energy matters.

## Trace a native pressure calculation

Consider a user program that constructs `PengRobinson76` and calls
`model%pressure(n, V, T, P)`. Follow this sequence in the source:

1. **Construction:** `src/models/residual_helmholtz/cubic/implementations/implementations.f90`
   defines `PengRobinson76`. It configures a `CubicEoS` using component data,
   pure-component parameters, an alpha function, and a mixing rule.
2. **Shared property:** `pressure` in
   `src/models/residual_helmholtz/ar_models.f90` requests the volume derivative
   of residual Helmholtz energy and forms the pressure relation.
3. **Dynamic dispatch:** the call to the model's `residual_helmholtz` binding
   reaches the cubic implementation in `cubic/generic_cubic.f90`. The concrete
   type, not the declared `class(ArModel)` variable, determines that binding.
4. **Model composition:** the cubic implementation evaluates its alpha and
   mixture parameters, then calls the algebraic energy kernel in
   `cubic/base.f90` for the requested quantities.

The same user pressure call on a hyperdual model goes through its AD adapter
instead of the cubic energy kernel. The shared pressure identity does not
change. This is the practical meaning of polymorphism in this library.

## Trace an equilibrium calculation

For a saturation pressure, begin in `src/equilibria/saturations_points.f90`
(note the existing filename spelling). The routine organizes initialization,
phase/root selection, iterations, and construction of `EquilibriumState`.
Lower-level simultaneous saturation equations live in
`src/models/solvers/saturation_point.f90`.

Those equations call model property methods. They do not contain a separate
Peng–Robinson thermodynamic implementation. If a residual is wrong, first
identify whether the error is in unpacking the solver variables, calling the
right phase/root, evaluating model derivatives, or assembling the residual.
Changing a generic Newton tolerance is rarely the first useful diagnosis.

Phase boundaries add continuation: predict a nearby state, correct the
equilibrium equations, choose a specification, and store the point. The
`src/equilibria/boundaries/` modules contain several specialized formulations;
not all of them use the same driver or unknown-vector layout. Read the
specific routine rather than assuming a single universal envelope algorithm.

## The language boundary is another layer

Native polymorphic model objects do not pass directly into Python. The
`yaeos_c` module stores model objects behind integer handles and offers flat
procedures. f2py generates a Python extension around that module; Python
classes hold the handle and present a more idiomatic user API.

The native `ArModel` class and the Python `ArModel` class have related purposes,
but they are not the same object in two syntaxes. Shapes, option flags, result
containers, and memory lifetime need explicit adaptation. Read the
[C-facing layer](c_interface.html) and [Python layer](python_interface.html)
before changing a wrapper signature.

## Choose the right place for a change

| Proposed change | Likely starting point | Check downstream |
|-----------------|-----------------------|------------------|
| New cubic alpha correlation | `cubic/alphas/` | Constructor/configuration, derivatives, wrappers |
| New mixture combining rule | `cubic/mixing_rules/` | Composition/state derivatives and pure limits |
| Entirely new EOS energy expression | Model implementation and derivative adapter | Shared property use, volume initialization, exports |
| Incorrect flash material balance | Equilibrium residual/phase reconstruction | Result conventions and cross-language wrappers |
| Wrong NumPy output shape | Python method and flat wrapper | Fortran dimensions and memory ordering |
| Missing argument explanation | User guide and defining procedure's FORD comments | Generated API and Python documentation |

## A short source-reading exercise

Choose a property you already use. Find its public export, defining procedure,
the model method it calls, and a regression test. Write down the units and
shape at each step. Then repeat through Python and mark where its arguments
are converted. This small exercise gives a more useful mental model than
reading the repository alphabetically.

Continue with [Fortran contracts](fortran.html), or use the
[end-to-end change walkthrough](walkthrough.html) to turn that map into a
contribution plan.
