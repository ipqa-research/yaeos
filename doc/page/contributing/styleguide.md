---
title: Style guide
---

# Readability is part of the numerical interface

Thermodynamic code is easier to review when readers can connect each quantity
to an equation, its units, and its independent variables. Prefer descriptive
names for algorithms and conventional symbols for equations; define the latter
where they first appear. Use three-space indentation, `implicit none`, explicit
imports where practical, and named closing statements.

## One variable per declaration

In new code and pedagogical examples, put **one variable on each declaration
line**, including scalars, arrays, loop indices, and derived-type components.
This leaves room for a useful unit or physical meaning and makes changes easy
to review. For example:

```fortran
use yaeos, only: pr
implicit none

real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
real(pr), intent(in) :: V    !! Total volume [L], not molar volume
real(pr), intent(in) :: T    !! Temperature [K]
real(pr), intent(out), optional :: ArV !! (dAr/dV) at fixed n,T [bar]
real(pr), allocatable :: kij(:, :) !! Dimensionless interaction matrix
integer :: i !! Component index
integer :: j !! Component index
```

This is a declaration fragment inside a procedure, not a standalone program.
Do not group these as `real(pr) :: n(2), V, T` or `integer :: i, j`. The
repository guide included below contains older grouped-array examples; use the
one-variable convention for new declarations rather than copying those groups.
Do not reformat unrelated source as part of a small fix.

The same rule applies to derived-type components and to dummy arguments:

```fortran
type :: MyModel
   real(pr), allocatable :: ac(:) !! Pure-component attraction parameters [bar L^2/mol^2]
   real(pr), allocatable :: b(:)  !! Pure-component covolumes [L/mol]
contains
   procedure :: residual_helmholtz => my_ar !! Evaluate A^r and its derivatives
   procedure :: get_v0 => my_v0 !! Volume scale for the liquid root [L]
end type MyModel
```

A reader comparing your code with the [Fortran core guide of the developer
handbook](../developer/fortran.html) should find the same conventions in
both.

Comments should add information: “total volume [L]” is useful; “variable V” is
not. Array comments should explain index order when it matters, especially for
interaction matrices and Hessians.

## Precision, units, and derivative names

- Import `pr` and use `real(pr)` and kind-qualified real literals such as
  `1.0_pr`. Import `R` instead of hard-coding a gas constant with different units.
- Use the library's bar, L, K, mol basis. An extensive energy is in bar L, not
  bar L/mol. State any conversions at external data boundaries.
- Use thermodynamic input order `n, V, T` or `n, P, T`. Name optional outputs
  at call sites when positional lists would obscure the derivative requested.
- Energy derivatives use `ArV`, `ArT`, `ArTV`, `ArV2`, `Arn`, `Arn2`, etc.;
  general properties use names such as `dPdV`. State what stays fixed rather
  than relying on a short variable name alone.
- The shape of a mole gradient is `size(n)`; the mole Hessian has two component
  indices. A normalized composition derivative is not a mole derivative.

## Match contracts rather than improving them accidentally

Specify `intent` for procedure arguments when defining a new interface. When
overriding an existing deferred binding, however, match its declared attributes
exactly. The hyperdual `Ar` passed-object dummy has no `intent`; adding one can
break the override. Analytic `residual_helmholtz`, hyperdual `Ar`, and Tapenade
`ar` are different contracts, despite related names.

Guard optional outputs with `present`. Calculate shared working quantities in
local variables, not in an optional output that a caller may omit. Every
requested output must be assigned regardless of which others are present.
Document required initialization, array sizes, admissible domains, and any
assumptions such as symmetric parameter matrices.

## Explain numerical choices and test them

Explain why a bound, step size, or tolerance is appropriate and give its units
or scale. Avoid unexplained “magic” constants. Test ordinary states before
edge cases, then add dilute, pure-component, and near-boundary cases relevant
to the model. Do not turn a failed convergence check into a successful result
by silently replacing nonfinite values.

For documentation, pair a mathematical definition with one minimal call and
an interpretation of its output. Preserve existing examples and references
when expanding a page. Compile runnable snippets against the current interface
and label sketches honestly. See the [contributing workflow](index.html) for the
build, test, and FORD preview steps, and the
[developer handbook](../developer/index.html) for the surrounding codebase.

## Link to the rendered documentation

In guide text, link to the generated page, for example
`[Mixing rules](../usage/eos/cubics/mixing.html)`, rather than a local `.md`
file. Resolve the path relative to the rendered page, and check any `#anchor`
against the generated HTML. FORD does not automatically rewrite ordinary
Markdown links from `.md` to `.html`.

This differs from FORD's `ordered_subpage` front matter and include directives:
those identify **source files** and must retain their `.md` filenames.
An explicitly labelled external GitHub source-file link can also legitimately
end in `.md`; it is not an internal documentation navigation link.

## Repository-wide style reference

{!STYLE_GUIDE.md !}
