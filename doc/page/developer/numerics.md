---
title: Working on numerical algorithms
---

[TOC]

# From thermodynamic equations to a reliable solver

A phase calculation combines several independently fallible pieces: a physical
model, a choice of unknowns, residual equations, derivatives, a nonlinear
algorithm, and result interpretation. Treat these as separate layers when
debugging. Otherwise, a unit error can look like a convergence problem, and a
small residual can hide convergence to the wrong branch.

This page is a source-reading and modification guide. The
[user equilibrium pages](../usage/phase_equilibrium/index.html) explain which
calculation to choose; the [testing guide](testing.html) explains how to verify
an implementation.

## 1. Write the problem down before editing the loop

For the routine you are changing, identify:

* **Fixed inputs:** feed, pressure, temperature, or other constraints.
* **Unknown vector:** exactly which physical quantity each entry represents.
* **Transformations:** for example log pressure instead of pressure, or a
  phase-composition ratio instead of a mole fraction.
* **Residuals:** one equation per unknown, including normalization and
  continuation constraints.
* **Acceptance rules:** residual tolerance, valid domain, iteration limit,
  branch tests, and stability checks.
* **Returned state:** which arrays and flags are meaningful on failure.

Search for both the routine that packs `X` and the routine that unpacks it.
Index positions are part of the internal contract. PT, PX, TX, three-phase,
and general multiphase equations do not all use the same layout.

### Log coordinates change the Jacobian

Suppose \(u=\ln P\). A residual derivative with respect to the packed unknown
is not simply its pressure derivative:

\[
\frac{\partial F}{\partial u}=P\frac{\partial F}{\partial P}.
\]

Log coordinates keep pressure positive, but they do not by themselves make the
iteration well scaled or prevent overflow. A finite-difference check must
perturb `u`, not physical pressure, if the Jacobian is defined with respect to
the packed vector.

## 2. Know which numerical implementation you are calling

Several solver families coexist. Similar routine names do not imply identical
interfaces or guarantees.

| Source | Role to inspect |
|--------|-----------------|
| `src/math/newton.f90` | Configurable Newton solver with diagnostics and safeguarding |
| `src/math/nonlineal.f90` | Another multivariate Newton implementation and homotopy helpers |
| `src/math/math.f90` | Scalar Newton, numerical helpers, and MINPACK adapters |
| `src/math/linalg.f90` | Dense linear-system and eigenproblem wrappers |
| `src/math/continuation.f90` | Predictor/corrector driver with callback interfaces |
| `src/models/solvers/saturation_point.f90` | Saturation equations and local point correctors |
| `src/equilibria/` | Problem-specific loops and phase logic |

When reading an unqualified `call newton(...)`, follow its `use` association
to the actual module. Before replacing a solver, compare callback signatures,
output status conventions, scaling, maximum-iteration behavior, and which
quantities the caller assumes were filled on return.

Do not infer that a low-level linear solve propagates all LAPACK errors simply
because it calls LAPACK. Check whether the wrapper exposes or handles the
status. If you need stronger error reporting, plan that as an explicit
interface change with tests rather than quietly changing a caller's assumed
semantics.

## 3. Verify the Jacobian independently

For \(F(X)=0\), Newton's method solves
\(J\Delta X=-F\), where \(J_{ij}=\partial F_i/\partial X_j\). A useful
central-difference check of column \(j\) is

\[
J_{:,j}\approx\frac{F(X+h_j e_j)-F(X-h_j e_j)}{2h_j}.
\]

Use a relative and absolute scale for \(h_j\) so a zero coordinate still gets
a nonzero perturbation. Try more than one step size: an extremely small step
can amplify roundoff, while a large step measures curvature rather than the
local derivative. Avoid crossing invalid domains or switching volume branches
when comparing derivatives intended for one smooth branch.

If the Jacobian is wrong, inspect unpacking, chain-rule factors, component
ordering, and phase conventions before changing damping. In particular, the
user convention \(K_i=y_i/x_i\) may be represented internally as an incipient
phase divided by the feed phase; dew calculations can therefore use inverted
ratios. Check the actual reconstruction equations.

## 4. Understand continuation as a sequence of local solves

A phase boundary is a curve, not a single root. A continuation method fixes
one specification \(S\) to select a point, corrects the equilibrium equations,
then predicts a nearby point using a tangent. For equations \(F(X,S)=0\),

\[
J\frac{dX}{dS}=-\frac{\partial F}{\partial S}.
\]

The generic driver in `src/math/continuation.f90` exposes callbacks for
residuals/derivatives, optional specification updates, postprocessing, point
correction, and stopping. Read the call site to see which callbacks a given
envelope actually supplies.

The current generic driver's default specification update selects the largest
absolute tangent component and rescales the tangent and step. It does **not**
apply the `0.05` step cap shown in an example comment. If a branch needs a
physical step limit, identify the active callback that supplies it or add a
tested policy deliberately.

The generic result has a fixed row count, with unfilled rows initialized to
zero. That is a storage convention, not proof that every requested point was
accepted. Specialized envelope routines can use other storage conventions.
Do not interpret padding as a zero-temperature equilibrium point.

## 5. Separate three meanings of “good result”

1. **Numerical convergence:** the algorithm met its stated termination test.
2. **Physical admissibility:** compositions, phase fractions, volumes, and
   temperatures lie in the intended domain and balances close.
3. **Physical relevance:** the branch is stable or otherwise appropriate for
   the task, and the model/parameters are valid for the system.

A converged negative phase fraction can be useful in a mathematical
continuation formulation but is not an ordinary physical phase split.
Likewise, forcing a liquid root may intentionally follow a metastable branch.
Preserve the distinction in names, documentation, tests, and wrappers.

Near a critical point, phases become indistinguishable and the equations may
be nearly singular. Tolerances and step control matter, but singularity can
also reflect the formulation itself. A full mixture criticality condition
includes the ideal contribution to the stability matrix; a residual-energy
Hessian alone is not the general mixture condition.

## 6. Build a minimal failure case

When an envelope fails at point 137, save the last accepted state, the
predictor, selected specification, residual norms, and relevant parameters.
Try to reproduce the **single point correction** before rerunning a long
curve. This makes it easier to tell whether the problem began with a bad
predictor, an invalid phase reconstruction, or a broken residual/Jacobian.

For a flash regression, check the feed balance in addition to a solver flag:

\[
z_i=(1-\beta)x_i+\beta y_i.
\]

For multiple phases, use the representation and reference-phase convention of
that routine rather than extrapolating the two-phase array layout. Never
compare two implementations' raw packed vectors until you have checked that
they pack the same variables in the same order.

## 7. Review numerical changes as behavioral changes

Changing a tolerance, line search, initial guess, or root choice can change
which solution is found. Record why the change is needed and test more than
the failing example: one ordinary case, the regression, a relevant limiting
case, and a failure path. For performance changes, compare both outputs and
timings under the same compiler/profile.

Avoid broadly increasing iteration limits to hide unexplained regressions.
Preserve a minimal test of the failure and document any remaining limitation.
The [end-to-end walkthrough](walkthrough.html) shows how to carry that evidence
through a contribution.
