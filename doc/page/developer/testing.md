---
title: Tests and numerical validation
---

[TOC]

# Write a test that catches the bug you fixed

A useful test distinguishes correct behavior from a plausible wrong answer.
For `yaeos`, that means testing mathematical contracts, units, array shapes,
phase balances, and failure handling—not only checking that an iteration
finished or a value is finite.

This guide starts with a small native test, then adds the independent checks
needed for thermodynamic and cross-language work. Use
[workflow](workflow.html) for prerequisites and build commands.

## Where tests live and how they run

The native tests are under `test/`; fixtures and model-specific cases appear
in its subdirectories. **Write new tests as standalone Fortran programs.**
The project is migrating away from test-drive suite registration. Older suite
modules and `test/test_runner.f90` remain in the tree during that migration,
but they are not the pattern to copy for a new test.

`fpm` discovers test programs automatically because `auto-tests = true` in
`fpm.toml`. Each program is a separate test target. Shared helpers and fixtures
can remain modules, but a module alone is not an executable test. Use
`fpm test --list` to confirm that your new program is discovered.

Python tests are under `python/tests/`. Some compare directly with native
Fortran calls by generating small programs in `app/`. Those tests need a
checkout, `fpm`, and a working native build as well as the Python package.

## Add a small native regression

### 1. Write a standalone program

Save this complete program as `test/test_my_feature.f90` in a working
development branch. It follows the program-based pattern used by
`test/test_huron_vidal.f90`, `test/test_cmr.f90`, and `test/test_sddlc.f90`:

```fortran
program test_my_feature
   use yaeos, only: pr, ArModel, PengRobinson76
   use testing_aux, only: assert, test_title
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   class(ArModel), allocatable :: model
   real(pr) :: n(2)
   real(pr) :: specified_pressure
   real(pr) :: calculated_pressure
   real(pr) :: volume
   real(pr) :: temperature

   write(*, *) test_title("Volume/pressure regression")
   model = PengRobinson76([190.564_pr, 305.32_pr], &
      [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr])
   n = [0.3_pr, 0.7_pr]
   specified_pressure = 10.0_pr
   temperature = 320.0_pr

   call model%volume(n, specified_pressure, temperature, volume, &
      root_type="stable")
   call assert(ieee_is_finite(volume), "Volume must be finite")
   call assert(volume > 0.0_pr, "Volume must be positive")

   call model%pressure(n, volume, temperature, calculated_pressure)
   call assert(ieee_is_finite(calculated_pressure), "Pressure must be finite")
   call assert(abs(calculated_pressure - specified_pressure) < 1.0e-7_pr, &
      "Volume reproduces the specified pressure")
end program test_my_feature
```

`testing_aux` is the repository's test helper module, defined in
`test/yaeos_testing.f90`, not part of the installed public `yaeos` API.
`test_title` formats a heading; `assert(condition, message)` reports a passing
check or terminates the program with `error stop` on failure. A failing
prerequisite therefore stops before later calculations use an invalid result.
Printing a warning without a nonzero exit status is not a failing test.

For a group of related checks, put private test subroutines after `contains`
and call them explicitly from the program body. There is no collector or
suite registry to update.

### 2. Let fpm discover and run the program

From the repository root:

```bash
fpm test --list
fpm test --target test_my_feature
```

The target comes from the source filename, so use a descriptive `test_*.f90`
name and a matching program name. No change to `test_runner.f90` or the
manifest is needed for an automatically discovered test program. Run the
focused target first, then `fpm test` before review. If a target is missing,
check that its source contains a `program` unit rather than only a module.

## Decide what kind of evidence a test provides

| Check | What it establishes | What it does not establish |
|-------|----------------------|----------------------------|
| Compare with independent data | Agreement for specified states and parameters | Accuracy everywhere |
| Numerical derivative check | Local agreement of derivatives and energy evaluation | That the energy equation is the intended physical model |
| Thermodynamic identity | Internal consistency of related quantities | Independence of shared implementation errors |
| Material balance and isofugacity | A candidate phase split satisfies its equations | Global stability against all other splits |
| Native/Python comparison | Agreement across wrapper paths | Independent verification of their shared native kernel |
| Regression reproducer | A particular failure remains fixed | Coverage of unrelated edge cases |

Use more than one kind of evidence for a new model or solver. Two routines
that share the same incorrect derivative can agree perfectly.

## Derivatives: test the variables the interface actually uses

For a scalar property \(f(X)\), a central-difference estimate is

\[
\frac{\partial f}{\partial X_j}\approx
\frac{f(X+h_j e_j)-f(X-h_j e_j)}{2h_j}.
\]

Choose steps with both an absolute and relative scale so that a zero input
still gets a nonzero perturbation. Repeat at several step sizes. Small steps
amplify cancellation; large steps introduce truncation error. Neither choosing
the smallest possible number nor copying another test's step blindly is a
reliable validation strategy.

Test the **amount derivatives** required by `ArModel` and `GeModel`, not
derivatives of normalized mole fractions. For an equilibrium Jacobian in log
coordinates, perturb the packed log variable and include its chain-rule
factor. Avoid changing volume-root branches during a derivative check intended
for a single smooth branch.

Useful native helpers include `numeric_ar_derivatives`, `ar_consistency`,
`numeric_ge_derivatives`, `ge_consistency`, and `individual_ar_calls` under
`src/consistency/`. Read their interfaces and the cases that call them:
different outputs have different scalar, vector, and matrix shapes.
`individual_ar_calls` helps check optional-output combinations; it does not
replace a mathematical accuracy check.

For a derivative implementation, cover energy only, each derivative separately,
all outputs together, mixed-output requests, and repeated calls at different
states. AD implementations need these checks too: differentiation propagates
the program's operations, including mistakes in its primal expression.

## Thermodynamic checks that expose amount or reference errors

### Scaling and homogeneity

For residual Helmholtz energy at fixed temperature,

\[
A^r(\lambda\mathbf n,\lambda V,T)=\lambda A^r(\mathbf n,V,T).
\]

Pressure and fugacity coefficients should stay unchanged under that same
joint scaling. Scaling amounts without volume changes the density, so it is
not the same test.

For an excess-Gibbs model,
\(G^E(\lambda\mathbf n,T)=\lambda G^E(\mathbf n,T)\), whereas activity
coefficients are unchanged. These checks help catch a missing total-amount
factor or an inappropriate normalization.

### Heat-capacity identity

At fixed component amounts, define
\(P_T=(\partial P/\partial T)_{V,\mathbf n}\) and
\(P_V=(\partial P/\partial V)_{T,\mathbf n}\). The total heat capacities
satisfy \(C_p-C_v=-T(P_T)^2/P_V\). Subtracting the ideal-gas difference gives

\[
C_p^r-C_v^r=-T\frac{(P_T)^2}{P_V}-n_{\mathrm{tot}}R.
\]

Here \(P_T=n_{\mathrm{tot}}R/V-A^r_{VT}\) and
\(P_V=-n_{\mathrm{tot}}RT/V^2-A^r_{VV}\). First derivatives of residual
Helmholtz energy cannot simply be multiplied to obtain this identity.

This **insertable test fragment** assumes an initialized model and declared
`n`, `V`, and `T`; import `pr` and `R` from `yaeos` and `assert` from
`testing_aux`. Place the declarations before executable statements.
Choose a regular state where `dPdV` is nonzero:

```fortran
real(pr) :: P
real(pr) :: dPdT
real(pr) :: dPdV
real(pr) :: Cp
real(pr) :: Cv
real(pr) :: expected_difference

call model%pressure(n, V, T, P, dPdT=dPdT, dPdV=dPdV)
call model%Cp_residual_vt(n, V, T, Cp=Cp)
call model%Cv_residual_vt(n, V, T, Cv=Cv)
expected_difference = -T * dPdT**2 / dPdV - sum(n) * R
call assert(abs(Cp - Cv - expected_difference) < 1.0e-8_pr, &
   "Residual heat capacities satisfy the pressure-derivative identity")
```

This tests consistency of property relations and derivative paths. An
independent derivative check is still needed to diagnose a shared model error.
Near singular states, conditioning and cancellation can make this a poor
choice of numerical test even when the identity is valid.

### Symmetry, Gibbs–Duhem, and limits

For a smooth excess-Gibbs model, the Hessian of \(G^E\) is symmetric and
\(\sum_i n_i\partial\ln\gamma_i/\partial n_j=0\). `ge_consistency` checks
these relations and the link between energy and activity coefficients. The
matrix of first mole derivatives of \(\ln\gamma\) is a **Jacobian**, not a
Hessian of \(\ln\gamma\).

Check dilute-gas limits where applicable: \(Z\to1\), \(\ln\phi_i\to0\),
and residual energy tends to zero as density tends to zero at fixed amount
and temperature. For activity models using pure-component reference states,
check appropriate pure limits and parameter conventions. Use a sufficiently
dilute finite state, not a singular zero-volume or zero-amount calculation.

## Phase calculations: check physics as well as termination

For an accepted two-phase result, check:

* Finite, nonnegative phase compositions, each summing to one.
* Finite phase fractions in their intended domain.
* Component balance: `z = (1-beta)*x + beta*y`.
* Volume balance for a VT flash, on the routine's specified molar basis.
* Log-fugacity agreement on the same roots used by the solver.
* Distinct phases when the requested calculation should not be trivial.

For multiphase output, reconstruct the feed using **all** phase fractions,
including the reference phase. `MPEquilibriumState%np` counts main phases in
the current representation; read the corresponding routine's contract rather
than assuming it is the total count.

Test ordinary cases, relevant boundaries, and a failure path. Keep
near-critical conditioning tests distinct from ordinary accuracy tests. A
converged root can be metastable, and a finite search that finds no destabilizing
trial phase is not a proof of global stability.

Failure conventions differ. Some routines return a sentinel, some return an
early result or status field, and some terminate with `error stop`. Test the
specific contract; do not read uninitialized result arrays after a failed
call. To test an **expected** terminating failure, use a subprocess and check
its exit status from a harness. Do not add a program that always exits nonzero
as an ordinary passing `fpm test` target: `fpm` correctly treats that exit as
a test failure.

## Cross-language checks: shapes before values

Python tests under `python/tests/` include native-reference comparisons. Some
generate Fortran programs in `app/`, execute them with `fpm`, and parse their
output with helpers such as `parse_vec` and `parse_mat`. Inspect an existing
test for the exact protocol instead of inventing another incompatible one.

For an interface change, check:

1. The same model parameters, state, units, and component ordering on both sides.
2. Scalar/vector/matrix shape before comparing values.
3. Matrix orientation with a nonsymmetric test case where appropriate.
4. Dtype, contiguous layout, and Fortran-order buffers where the wrapper requires them.
5. Each optional derivative alone and in combination.
6. The return container and unrequested-output convention, not only its numbers.
7. Model creation/deletion and invalid-handle behavior where safely defined.

A Python test that agrees with native Fortran validates the adaptation layer,
not the physical model independently. Test an installed build after changes to
the compiled wrapper; a stale extension can conceal or mimic a regression.

## Choose tolerances to match the claim

A practical floating-point comparison uses both scales:

\[
|a-b|\leq \mathrm{atol}+\mathrm{rtol}|b|.
\]

Choose `atol` in the quantity's units and `rtol` for its expected scale.
Near-zero consistency residuals need an absolute criterion. Properties spanning
orders of magnitude often need a relative criterion too. Do not compare
log-fugacity mismatch, temperature error, and material balance with a single
unexplained tolerance.

State the reason for a reference value and its tolerance in the test. Consider
rounding in published parameters, conditioning, compiler differences, and
finite-difference error. Do not loosen a tolerance until a failure disappears
without understanding why.

## A minimal pre-review test plan

Run the affected target first, then the native suite. Rebuild the Python
extension and run relevant Python tests for an interface change. Build the
documentation and run or compile its complete examples separately: successful
site generation does not execute the tutorial.

In the review summary, report which layers were tested, the revision and
environment, and anything unavailable. Coverage is useful for finding unvisited
branches, but it cannot establish that the equations, reference data, or
assertions are correct.

Continue with [workflow](workflow.html), [numerical implementation](numerics.html),
or [the end-to-end change walkthrough](walkthrough.html).
