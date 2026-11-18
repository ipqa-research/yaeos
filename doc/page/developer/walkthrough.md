---
title: An end-to-end change
---

[TOC]

# Plan a small feature from calculation to public API

This walkthrough uses a small property helper to show how a contribution
grows across the codebase. The mathematical task is deliberately simple:
calculate the compressibility factor

\[
Z=\frac{PV}{n_{\mathrm{tot}}RT}.
\]

\(Z\) is dimensionless and equals one for an ideal gas. It is not the same
quantity as the isothermal compressibility, which is a derivative of volume
with respect to pressure. Naming the desired property precisely is the first
part of a good feature proposal.

The helper below is **application code**, not a new method already supplied by
`ArModel`. No library modification is necessary to run it. The later sections
describe a possible contribution plan, not commands that silently add that
method to your checkout.

## 1. Define a contract in a few sentences

Given an EOS model, component mole amounts `n`, total volume `V` in L, and
temperature `T` in K, evaluate model pressure in bar and return \(Z\). Use the
total amount, not the length of the composition array or an implicit one-mole
normalization. The physical state is homogeneous: the helper does not perform
a flash or decide whether the mixture is stable.

Before adding a method to the library, search for equivalent functionality and
discuss whether a helper is useful enough to publish. A one-line application
calculation can be preferable to expanding a maintained public interface.

## 2. Make the calculation work outside the library first

The following is a complete Fortran source file. Put it in an application's
`app/main.f90` after configuring the `yaeos` dependency as in
[getting started](../usage/index.html):

```fortran
module application_properties
   use yaeos, only: pr, R, ArModel
   implicit none
   private
   public :: z_factor
contains
   function z_factor(model, n, V, T) result(Z)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: V
      real(pr), intent(in) :: T
      real(pr) :: Z
      real(pr) :: P

      if (any(n < 0.0_pr) .or. sum(n) <= 0.0_pr) &
         error stop "Invalid mole amounts"
      if (V <= 0.0_pr .or. T <= 0.0_pr) error stop "Invalid V or T"
      call model%pressure(n, V, T, P)
      Z = P * V / (sum(n) * R * T)
   end function z_factor
end module application_properties

program check_property_helper
   use yaeos, only: pr, ArModel, PengRobinson76
   use application_properties, only: z_factor
   implicit none

   class(ArModel), allocatable :: model
   real(pr) :: n(2)
   real(pr) :: Z
   real(pr) :: Z_scaled
   real(pr) :: Z_dilute

   model = PengRobinson76([190.564_pr, 305.32_pr], &
      [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr])
   n = [0.3_pr, 0.7_pr]
   Z = z_factor(model, n, 1.0_pr, 320.0_pr)
   Z_scaled = z_factor(model, 10.0_pr * n, 10.0_pr, 320.0_pr)
   Z_dilute = z_factor(model, n, 10000.0_pr, 320.0_pr)

   if (abs(Z - Z_scaled) > 1.0e-10_pr) error stop "Scaling check failed"
   if (abs(Z_dilute - 1.0_pr) > 1.0e-3_pr) &
      error stop "Dilute-gas check failed"
   write(*, '(A,F12.8)') "Z at the test state: ", Z
   write(*, '(A,F12.8)') "Z in the dilute limit: ", Z_dilute
end program check_property_helper
```

The `error stop` checks are convenient in this demonstration application. They
are not a recommendation to terminate every host process from a public
library routine. A production API needs deliberate invalid-input and failure
reporting, especially when called from Python or another language.

## 3. Choose tests that can catch your own mistakes

This example checks two independent properties:

* **Scaling:** multiplying amounts and total volume by ten must leave \(Z\)
  unchanged. This catches accidental normalization of only one input.
* **Dilute behavior:** at a sufficiently large volume, a suitable EOS should
  approach ideal-gas behavior. This catches some unit and formula errors.

Neither check proves the EOS is accurate for methane–ethane at all conditions.
Add an ordinary regression state and a known reference calculation if the
property becomes public. Also test the one-component case, different amount
bases, component reordering with corresponding model reordering, and relevant
invalid inputs. Use an appropriate tolerance for each claim, not one absolute
tolerance for every physical quantity.

For a bug fix, preserve the smallest input that previously failed. Keep the
regression test after the fix, and explain the failure mechanism in the change
description. See [testing](testing.html) for native test-program and Python
test conventions.

## 4. Decide where the core feature belongs

There are two plausible designs: a standalone utility that accepts
`class(ArModel)`, or a shared type-bound method that uses the existing pressure
binding. Neither needs a separate implementation for each cubic EOS.

If choosing a type-bound method, inspect
`src/models/residual_helmholtz/ar_models.f90` for the surrounding conventions.
If choosing a standalone public utility, check the relevant facade exports.
Document units, the amount basis, whether pressure is computed or supplied,
and whether a volume-root choice exists.

Do not simultaneously change all neighboring code to a new style. Keep the
functional change small enough that a reviewer can distinguish its contract
from unrelated refactoring.

## 5. Add language exposure only when it is needed

A Python-only convenience method can already calculate this quantity using
the existing pressure wrapper:

```python
import numpy as np
from yaeos import PengRobinson76

model = PengRobinson76(
    [190.564, 305.32], [45.99, 48.72], [0.01142, 0.0995]
)
moles = np.array([0.3, 0.7])
volume = 1.0
temperature = 320.0
pressure = model.pressure(moles, volume, temperature)
gas_constant = 0.08314462618  # bar L / (mol K), matching the Fortran library
z = pressure * volume / (moles.sum() * gas_constant * temperature)
```

This is a complete Python example when an appropriate `yaeos` installation is
available. A user-level convenience calculation does not automatically justify
another compiled entry point.

For a feature that **does** need a new native wrapper, make an explicit plan:

| Layer | Change to consider | Evidence to add |
|-------|--------------------|-----------------|
| Fortran kernel/shared property | Implement the mathematical contract | Native regression and limiting cases |
| Fortran public API | Add binding/export if needed | A small external caller compiles |
| `yaeos_c` wrapper | Adapt handles, kinds, dimensions, and optional results | Wrapper matches native values and shapes |
| Python method | Convert inputs and preserve result/derivative conventions | Python regression and invalid-input tests |
| Package/build metadata | Export or install new modules if necessary | A fresh installation imports them |
| FORD and Sphinx guides | Explain usage, restrictions, and references | Both documentation builds and links checked |

Use existing wrappers as patterns, but verify their assumptions. Avoid adding
a supposedly C-callable routine solely by using `real(c_double)`: interoperability
also depends on binding, shapes, argument passing, and ownership. See the
[C-facing interface guide](c_interface.html).

## 6. Review compatibility before changing an existing signature

Adding a native optional argument can affect a generated wrapper's signature,
even if old Fortran calls still compile. Changing array orientation can
silently transpose a result in Python. Renaming a result field can break
plotting and fitting tools that never call the Fortran routine directly.

Search both native and Python call sites, examples, tests, and documentation.
For a new model type, also check constructors, facade exports, handle storage,
Python package exports, source-install lists in Meson, and object lifetime.

## 7. Submit evidence, not just a passing screenshot

Explain the intended behavior, the affected layers, the numerical basis of the
tests, and the validation commands actually run. State any unavailable check,
such as an unbuilt Python extension. A helpful review summary distinguishes
“native Fortran tests pass”, “Python wrapper tests pass”, and “documentation
examples compile”; none is a substitute for the others.

Use [workflow](workflow.html) for the build commands and
[the developer landing page](index.html) to choose the next area to learn.
