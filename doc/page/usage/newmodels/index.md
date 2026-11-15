---
title: Adding your own models
ordered_subpage: analtical.md
                 autodiff.md
---

Most thermodynamic properties calculated in `yaeos` depend heavily on different
kinds of models and their derivatives. Obtaining the derivatives of complex
models can be tedious and error-prone, so `yaeos` provides three ways to add a
new model:

1. Analytical derivatives, implemented by hand.
2. Automatic differentiation with hyperdual numbers.
3. Source-to-source differentiation with Tapenade.

These are implementation choices, not different thermodynamic models. Begin
with the same energy equation and choose how to obtain its derivatives.

[TOC]

# Configure first, extend only when necessary

You do **not** need a new derived type to change critical constants, interaction
parameters, or parameters of an existing activity model. Start with the
[equations of state](../eos/index.html) and
[excess Gibbs models](../excessmodels/index.html). For a cubic EOS, an existing
alpha function or mixing rule may already express your intended change.
Extending just that component is usually smaller and easier to test than
rewriting the whole EOS.

Implement a new model when the energy equation itself is missing, or when you
need a genuinely new functional form. Before doing so, you should be able to:

1. Construct an existing model and calculate pressure at specified `n, V, T`.
2. Distinguish total volume from molar volume and residual from ideal-gas energy.
3. Read a Fortran module, a derived type, a type-bound procedure, and an optional
   dummy argument. The [style guide](../../contributing/styleguide.html) explains
   the conventions used here.
4. State the new equation's parameters, units, validity range, and limiting cases.

A useful beginner path is to reproduce the
[analytic van der Waals example](analtical.html), test its derivatives, then
translate the primal energy into hyperdual arithmetic. Only then tackle a more
complex model or a Tapenade generation workflow. Keep the first tests at a
single, well-behaved state; phase-equilibrium solvers are not derivative debuggers.

# Why start from an energy and its derivatives?

An EOS provides an extensive residual Helmholtz energy
\(A^r(\mathbf n,V,T)\): the departure from the ideal gas at the same
`n, V, T`. For example,

\[
P=\frac{(\sum_i n_i)RT}{V}-A^r_V,\qquad S^r=-A^r_T.
\]

Mole derivatives supply residual chemical potentials; second derivatives supply
property sensitivities and the Jacobians used by nonlinear solvers. Deriving
these quantities from a common energy helps preserve thermodynamic identities.
Do not add the ideal-gas pressure inside `residual_helmholtz`: the inherited
`pressure` method already adds it.

An activity model instead supplies extensive excess Gibbs energy
\(G^E(\mathbf n,T)\), relative to an ideal solution, with
\(\ln\gamma_i=G^E_{n_i}/(RT)\). It is not an EOS and does not independently
provide pressure or a liquid density. Its abstract base is [[GeModel]];
implement its `excess_gibbs` subroutine with the optional energy, temperature,
and mole derivatives in `src/models/excess_gibbs/ge_models.f90`.
There is no `get_v0` obligation for a `GeModel`.

Both energies are returned in **bar L**, not joules, molar energy, or energy
divided by `RT`. Use `R` supplied by `yaeos`; 1 bar L is 100 J. Inputs are
component mole amounts [mol], total volume [L] for an EOS, and temperature [K].
Using mole fractions as `n` is a one-mole basis, not permission to normalize
`n` internally. A composition derivative changes one mole amount while holding
the other mole amounts fixed, not while forcing their sum to one.

# Residual Helmholtz models
Residual Helmholtz models are the basis to obtain the residual properties.

The main basis in `yaeos` is the abstract type [[ArModel]], which enforces the
expected functionality of this kind of model. **Abstract** means that it
describes a contract, not a model you can instantiate. A **concrete type** stores
parameters and supplies all deferred procedures in that contract. An abstract
interface specifies the procedure kind (function or subroutine), argument
names, types, ranks, and attributes. The compiler checks your bindings against
it; matching just the procedure name is not enough.

This deliberately incomplete sketch illustrates inheritance, not a compilable
concrete model:

```fortran
use yaeos, only: ArModel

type, abstract, extends(ArModel) :: MyNewModel
end type
```

The definition of an `ArModel` requires two procedures:

- [[abs_residual_helmholtz(interface)]]: calculates the residual Helmholtz
  energy and its derivatives.
- [[abs_volume_initializer(interface)]]: provides an initial volume to find a
  liquid root at a given pressure and temperature.

```fortran
use yaeos, only: ArModel

type, extends(ArModel) :: MyNewModel
contains
   procedure :: residual_helmholtz => an_Ar_implementation
   procedure :: get_v0 => an_v0_implementation
end type
```

The binding sketch above needs module procedures with the exact interfaces.
Once those procedures and a parameter constructor are defined, the model can
use inherited methods (the following is a usage fragment):

```fortran
use my_model, only: MyNewModel

type(MyNewModel) :: model

! Assuming the model parameters are set up
call model%pressure(n, V, T, P)
```

`type(MyNewModel)` declares that specific implementation;
`class(ArModel), allocatable` can hold any concrete extension of `ArModel`.
That polymorphism lets the same flash or property routine work with different
EOS implementations without knowing their parameter layout.

# A validation ladder

Before trusting a new model, check in this order:

1. **Primal equation:** units, parameter dimensions, admissible logarithms and
   denominators, pure-component and low-density limits.
2. **Derivatives:** finite differences over several step sizes, Hessian symmetry,
   and each optional output requested alone as well as together.
3. **Extensivity:** at fixed `T`, scaling both `n` and `V` by a positive factor
   must scale `Ar` by that factor and leave pressure unchanged. For `Ge`, scale
   `n` only: `Ge` scales and activity coefficients stay unchanged.
4. **Thermodynamic identities:** collect the consistency residuals rather than
   assuming that a call without outputs prints a pass/fail report.
5. **Solvers:** test volume roots before saturation, flash, or envelopes; a
   converged result does not establish that the physical equation is correct.

See the analytic guide for executable checks and how to interpret failures.

# Guides

- [Analytical derivatives](analtical.html) shows how to implement the
  `residual_helmholtz` and `get_v0` procedures by hand, with a complete
  van der Waals example.
- [Automatic differentiation](autodiff.html) shows how to obtain the derivatives
  with hyperdual numbers and with Tapenade.
- [Excess Gibbs models](../excessmodels/index.html) describes how to add an
  activity coefficient model.

# Where this fits in the codebase

These pages teach the model author's side of the work. The
[developer handbook](../../developer/index.html) covers the surrounding
codebase: the repository map, the Fortran contracts, the numerical algorithms,
the build and test workflow, and the language-interface layers. Read it when
you need to know *where* something lives or *how* to build and test the
project, not only how to write the energy equation.

A few handbook pages pair directly with this guide:

| Need | Handbook page |
|:--|:--|
| Where a model, solver, or interface lives | [Architecture](../../developer/architecture.html) |
| Exact interface rules and optional-output handling | [Fortran core](../../developer/fortran.html) |
| Turning a derivative check into a regression test | [Testing](../../developer/testing.html) |
| Build, test, and documentation commands | [Workflow](../../developer/workflow.html) |
| Debugging a solver that will not converge | [Numerical algorithms](../../developer/numerics.html) |

# Example models

The source tree contains complete implementations that can be used as
references:

- `src/models/residual_helmholtz/cubic/` for cubic equations of state.
- `src/models/residual_helmholtz/multifluid/gerg2008.f90` for a multifluid
  model with hyperdual automatic differentiation.
- `src/models/residual_helmholtz/saft/pc_saft.f90` for a Tapenade-generated
  model.
- `src/models/excess_gibbs/` for activity-coefficient models.

The [architecture page of the developer handbook](../../developer/architecture.html)
places these in the type hierarchy and explains how a user call reaches them.
The [testing handbook](../../developer/testing.html) shows how the shipped
models are checked, which is the pattern to copy for a new one.
