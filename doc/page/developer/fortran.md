---
title: Working in the Fortran core
---

[TOC]

# Understand the Fortran contracts before extending them

The most important part of adding a model is not choosing a filename. It is
satisfying the contract that shared property routines and solvers already
expect. A function that returns a plausible pressure at one state can still
have incorrect composition derivatives, an incompatible volume initializer,
or an energy normalized on the wrong amount basis.

Use this page alongside the [repository map](architecture.html). For complete
model implementations, continue to the
[analytical](../usage/newmodels/analtical.html) and
[automatic-differentiation](../usage/newmodels/autodiff.html) tutorials.

## Read the type hierarchy in ordinary Fortran terms

A **derived type** groups data and procedures. A type can extend a parent and
inherit its data and methods. A **deferred binding** specifies a procedure
interface but leaves its implementation to a concrete extension. A variable
declared with `class(ArModel)` can refer to different concrete EOS types.

```fortran
! Fragment inside a program that imports these names from yaeos.
class(ArModel), allocatable :: model
model = PengRobinson76(tc, pc, omega)
call model%pressure(n, V, T, P)
```

The assignment constructs a concrete model. The pressure routine then calls
the `residual_helmholtz` binding associated with that concrete type. You do not
need a `select type` branch for every EOS in a generic property calculation.
Add one only when the operation genuinely requires model-specific information.

`BaseModel` contains `components`, a `Substances` value whose arrays hold
critical properties and other metadata in component order. A constructor must
initialize the fields its algorithms actually use; leaving arrays unallocated
and hoping every inherited method tolerates that is not a valid contract.

## The residual Helmholtz interface

The defining interface is `abs_residual_helmholtz` in
`src/models/residual_helmholtz/ar_models.f90`. Its inputs are component mole
amounts `n(:)`, total volume `V`, and temperature `T`. Outputs are optional:

| Output | Meaning | Shape |
|--------|---------|-------|
| `Ar` | Residual Helmholtz energy | Scalar |
| `ArV`, `ArT` | First volume and temperature derivatives | Scalar |
| `ArV2`, `ArT2`, `ArTV` | Second state derivatives | Scalar |
| `Arn` | Mole-number gradient | `(nc)` |
| `ArVn`, `ArTn` | State/mole-number mixed derivatives | `(nc)` |
| `Arn2` | Mole-number Hessian | `(nc,nc)` |

Here `nc = size(n)`. The derivatives are with respect to **amounts**, not
independently varied mole fractions. Holding total volume fixed while varying
one amount changes both density and composition. Normalizing `n` before
differentiating can silently implement a different derivative.

### Extensive energy is a requirement

The energy is in bar·L, not automatically J/mol or the reduced molar quantity
\(A^r/(n_{\mathrm{tot}}RT)\). If a published correlation supplies a
dimensionless molar energy \(\alpha^r\), convert it consistently:

\[
A^r=n_{\mathrm{tot}}RT\alpha^r.
\]

Differentiate the full extensive expression, including its amount and
temperature factors. Adding those factors only to the final energy output
does not fix derivatives that were taken on a different function.

For fixed temperature, the extensive scaling identity is

\[
A^r(\lambda\mathbf n,\lambda V,T)=\lambda A^r(\mathbf n,V,T).
\]

It implies, for example,
\(A^r=VA^r_V+\sum_i n_i A^r_{n_i}\). These identities make useful tests
alongside independent derivative checks.

### Optional means every supported combination must work

Guard assignment to an optional output with `present`. Do not reference an
absent output as a temporary. This is a **schematic implementation fragment**,
not a complete energy kernel:

```fortran
! Compute local values independently of whether corresponding outputs exist.
if (present(ArV)) ArV = local_ArV
if (present(Arn)) Arn = local_Arn
```

Also consider what a caller requests internally. The current shared
`pressure` method requests several state derivatives even when the user asks
only for pressure. Conversely, direct calls can request just a Hessian or a
mixed derivative without requesting `Ar`. A model must not work only for the
combination used in its first example.

Tests should include energy only, individual derivatives, mixed combinations,
and all outputs together. Verify that repeated calls do not depend on
uninitialized work arrays or results left by a previous call.

## The excess-Gibbs interface

`GeModel` in `src/models/excess_gibbs/ge_models.f90` requires
`excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)`. The pattern is
similar, but there is no volume input and no EOS volume root in this contract.
`Gen` and `GeTn` are vectors; `Gen2` is a matrix.

The shared activity-coefficient method uses

\[
\ln\gamma_i=\frac{G^E_{n_i}}{RT}.
\]

As with Helmholtz energy, an excess-Gibbs implementation must return the
extensive energy on the library's unit basis. A molar correlation needs its
total-amount factor. Do not route a liquid activity model through an EOS
property calculation merely because both types inherit `BaseModel`.

## Choose a derivative implementation strategy

| Route | Implement | Good reason to choose it | Main risk to test |
|-------|-----------|--------------------------|-------------------|
| Analytical | The real-valued energy and derivative interface | Existing reliable formulas or optimized kernel | Algebra, index orientation, optional branches |
| Hyperdual | The adapter's hyperdual energy function and required model support | A clear primal expression with overloaded operations | Unsupported operations, dropped derivative lanes, nonsmooth branches |
| Tapenade | Primal routines plus generated derivative bindings | Source-to-source differentiation of a real-valued kernel | Seed/adjoint contracts, generated code integration and reproducibility |

Automatic differentiation applies chain rules to the program you wrote. It
does not know whether you implemented the intended thermodynamic equation or
used the right units. Keep analytical identities and independent numerical
checks in all three routes.

For an existing cubic model, first ask whether a new `AlphaFunction` or
`CubicMixRule` is enough. Replacing the entire `CubicEoS` energy implementation
to change one temperature correlation duplicates more logic than necessary.

## Volume initialization is a second contract

An `ArModel` also implements `get_v0(self, n, P, T)`. The default volume solver
uses it to initialize the dense branch. The specialized Michelsen volume
algorithm uses a covolume or lower-volume scale in a reduced-coordinate solve.
These roles are not identical.

Choose an initializer compatible with the actual volume method your model
uses. An ideal-gas volume returned for every branch is not automatically an
adequate liquid initializer. If a model overrides `volume`, check that the
override has the inherited signature and consistent root/failure behavior.

Read `src/models/solvers/volume.f90` and a comparable model before reusing that
solver. Test pressure recovered at the solved volume, extensive volume
scaling, and both requested root branches at appropriate states.

## Constructor, visibility, and publication

A new model module normally has three parts: a type that stores its parameters,
a constructor that initializes them, and the bound procedures. Keep parameter
validation near construction where possible. Document dimensions, units,
interaction matrix orientation, and whether inputs are copied.

To publish a model, trace the relevant imports through
`src/models/models.f90` and `src/yaeos.f90`. To publish an equilibrium routine,
also check `src/equilibria/equilibria.f90`. Use narrow imports inside
implementation modules to reduce unnecessary dependencies and name conflicts.

Adding a `use` line in the facade is not sufficient for Python support. The
flat wrapper, constructor handle management, Python class, and package exports
are separate changes. Follow the [end-to-end checklist](walkthrough.html).

## Write documentation that explains the contract

Use one variable declaration per line so each argument or work value can have
an unambiguous description:

```fortran
real(pr), intent(in) :: T !! Temperature [K].
real(pr), intent(out) :: P !! Pressure [bar].
```

For modules, types, and procedures, write a title and short summary, then
`# Description`, `# Examples`, and `# References` when references are needed.
Explain why the calculation exists and what a caller must guarantee, not just
the English expansion of its name. Mark schematic fragments as such; avoid
examples with undeclared variables presented as complete programs.

Source comments are part of the implementation's interface. If a change alters
units, failure behavior, shape, default settings, or derivative conventions,
update the comments, user guide, and relevant wrapper documentation together.

## A review checklist for a core change

Before review, be able to answer:

* Which binding or routine changed, and what callers depend on its signature?
* Does it preserve the amount basis, units, and derivative shapes?
* Are optional outputs and invalid inputs handled deliberately?
* Is the volume initializer compatible with the selected volume method?
* Are mathematical identities, individual derivatives, and regression cases tested?
* Are public exports and downstream interfaces updated if needed?

Continue with [numerical algorithms](numerics.html) if you are changing a
solver, or [testing](testing.html) to turn these contracts into regression
checks.
