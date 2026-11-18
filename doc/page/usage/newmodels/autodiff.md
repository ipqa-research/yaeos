---
title: Automatic differentiation
---

[TOC]

Implementing a new model and all its required derivatives can be tedious and
error-prone. Automatic differentiation (AD) avoids most of that work. Its
runtime cost depends on the equation, derivative requests, and implementation;
measure it before choosing an approach solely for speed.

Automatic differentiation applies the chain rule to the operations in a program,
so it computes derivatives of the implemented equations rather than
approximating them with finite differences. Forward mode propagates directional
derivatives from inputs to outputs; reverse mode propagates sensitivities back
from outputs to inputs. Combining modes gives mixed and second derivatives.
This is useful here because thermodynamic properties and phase-equilibrium
calculations require several derivatives of the residual Helmholtz energy.

In `yaeos` there are two ways to add a model with automatic differentiation:
operator overloading with `hyperdual` numbers and source-to-source
differentiation with Tapenade.

@warn
Remember to use the `R` constant from `yaeos__constants`. The base model already
provides a `type(Substances)` attribute named `components`; initialize it rather
than redeclaring it. Return extensive energy in bar L, not `Ar/(RT)` or molar
energy. Inputs are mole amounts [mol], total volume [L], and temperature [K].
@endwarn

## Hyperdual autodiff

### `ArModelAdiff`

Automatic differentiation with `hyperdual` numbers is available through the
[[ArModelAdiff]] derived type. A new model extends it and provides:

- `Ar`: a function returning the residual Helmholtz energy as a `hyperdual`
  value.
- `get_v0`: a volume initializer.

The interface of `Ar` is:

```fortran
type(hyperdual) function Ar(self, n, v, t)
   class(ArModelAdiff) :: self
   type(hyperdual), intent(in) :: n(:)  !! Number of moles
   type(hyperdual), intent(in) :: v      !! Volume [L]
   type(hyperdual), intent(in) :: t      !! Temperature [K]
end function Ar
```

This is an interface sketch. In a concrete implementation, replace
`class(ArModelAdiff)` by your concrete model type. Match the actual deferred
interface: its passed-object dummy has no `intent` attribute. In particular,
do not add `intent(in)` here just because the analytic API uses it.

The [[ArModelAdiff]] type implements `residual_helmholtz` by evaluating `Ar`
with hyperdual perturbations and extracting the first and second derivatives.
The model therefore gets all the derivatives required by the rest of the
library automatically.

### Step-by-step hyperdual workflow

1. Write down the extensive primal equation and test its units and domain with
   real arithmetic. AD differentiates what you wrote, not what you intended.
2. Extend `ArModelAdiff`, store real-valued parameters, and bind its `Ar`
   **function** and the inherited real-valued `get_v0` **function**.
3. Import `hyperdual_mod`, including its overloaded operators and intrinsics.
   Keep every intermediate depending on `n`, `v`, or `t` hyperdual-valued.
   Converting it to a real discards the derivative information. Model constants
   and inactive parameters may remain real.
4. Express the energy with supported operations such as `log`, `sqrt`, `sum`,
   and arithmetic. An external real-only helper will not automatically become
   differentiable: overload it or express its calculation in hyperdual form.
5. Let the adapter manage seeds; request properties through the ordinary
   `ArModel` API. Test each derivative and optional-output combination using
   the [validation ladder](analtical.html#validating-a-new-model).

Branches and clipping deserve particular care. A branch selects the derivative
of the executed path; a discontinuity does not acquire a meaningful smooth
derivative simply because AD is used. Test both sides of regime changes.

### Example

The following module implements a Peng-Robinson 1976 model with quadratic
mixing rules:

```fortran
module hyperdual_pr76
   use yaeos__constants, only: pr, R
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   use yaeos__substance, only: Substances
   use hyperdual_mod
   implicit none

   type, extends(ArModelAdiff) :: PR76
      !! Peng-Robinson 1976 model with hyperdual derivatives.
      real(pr), allocatable :: kij(:, :)
      real(pr), allocatable :: lij(:, :)
      real(pr), allocatable :: ac(:)
      real(pr), allocatable :: b(:)
      real(pr), allocatable :: k(:)
   contains
      procedure :: Ar => arfun
      procedure :: get_v0 => v0
   end type PR76

   real(pr), parameter :: del1 = 1._pr + sqrt(2._pr)
   real(pr), parameter :: del2 = 1._pr - sqrt(2._pr)

contains

   type(PR76) function setup(tc, pc, w, kij, lij) result(self)
      !! Set up the model from critical constants and interaction matrices.
      real(pr), intent(in) :: tc(:) !! Critical temperatures [K]
      real(pr), intent(in) :: pc(:) !! Critical pressures [bar]
      real(pr), intent(in) :: w(:)  !! Acentric factors
      real(pr), intent(in) :: kij(:, :) !! \(k_{ij}\) matrix
      real(pr), intent(in) :: lij(:, :) !! \(l_{ij}\) matrix

      self%components%tc = tc
      self%components%pc = pc
      self%components%w = w

      self%ac = 0.45723553_pr * R**2 * tc**2 / pc
      self%b = 0.07779607_pr * R * tc / pc
      self%k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2

      self%kij = kij
      self%lij = lij
   end function setup

   function arfun(self, n, v, t) result(ar)
      class(PR76) :: self
      type(hyperdual), intent(in) :: n(:)
      type(hyperdual), intent(in) :: v
      type(hyperdual), intent(in) :: t
      type(hyperdual) :: ar

      type(hyperdual) :: amix
      type(hyperdual) :: a(size(n))
      type(hyperdual) :: ai(size(n))
      type(hyperdual) :: nij
      type(hyperdual) :: bmix
      type(hyperdual) :: b_v
      integer :: i
      integer :: j

      associate( &
         ac => self%ac, b => self%b, k => self%k, &
         kij => self%kij, lij => self%lij, &
         tc => self%components%tc &
         )

         ! Soave alpha function
         a = 1.0_pr + k * (1.0_pr - sqrt(t/tc))
         a = ac * a**2
         ai = sqrt(a)

         ! Quadratic mixing rule
         amix = 0.0_pr
         bmix = 0.0_pr

         do j = 1, size(n) - 1
            do i = j + 1, size(n)
               nij = n(i) * n(j)
               amix = amix + 2 * nij * ai(i) * ai(j) * (1 - kij(i, j))
               bmix = bmix + nij * (b(i) + b(j)) * (1 - lij(i, j))
            end do
         end do

         amix = amix + sum(n**2 * a)
         bmix = bmix + sum(n**2 * b)
         bmix = bmix / sum(n)

         b_v = bmix / v

         ar = ( &
            -sum(n) * log(1.0_pr - b_v) &
            - amix / (R * t * bmix) * 1.0_pr / (del1 - del2) &
            * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
            ) * (R * t)
      end associate
   end function arfun

   function v0(self, n, p, t) result(v0out)
      !! Covolume initializer for the liquid root.
      class(PR76), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: v0out
      integer :: i
      integer :: j

      ! Same extensive covolume as in arfun, including lij.
      v0out = sum(n**2 * self%b)
      do j = 1, size(n) - 1
         do i = j + 1, size(n)
            v0out = v0out + n(i)*n(j)*(self%b(i) + self%b(j)) &
               * (1.0_pr - self%lij(i, j))
         end do
      end do
      v0out = v0out / sum(n)
   end function v0

end module hyperdual_pr76
```

Use symmetric interaction matrices with zero diagonals for this example and
check that the mixed covolume remains positive. The compact constructor assumes
matching array sizes and positive critical constants. `get_v0` returns the
same **total covolume** used by `arfun`, including `lij`. Dividing the linear
sum `sum(n*b)` by `sum(n)` would return a molar volume and only appear correct
on a one-mole basis. See the analytic guide for the default solver's `B/V`
transformation and domain obligations.

Once the model is defined, it can be used like any other `ArModel`:

```fortran
program main
   use yaeos, only: pr, ArModel
   use hyperdual_pr76, only: PR76, setup
   implicit none

   class(ArModel), allocatable :: model
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: kij(2, 2)
   real(pr) :: lij(2, 2)
   real(pr) :: n(2)
   real(pr) :: V
   real(pr) :: T
   real(pr) :: P

   tc = [190.0_pr, 310.0_pr]
   pc = [14.0_pr, 30.0_pr]
   w = [0.001_pr, 0.03_pr]

   kij = reshape([0.0_pr, 0.1_pr, 0.1_pr, 0.0_pr], [2, 2])
   lij = kij / 2.0_pr

   model = setup(tc, pc, w, kij, lij)

   n = [0.3_pr, 0.7_pr]
   V = 1.0_pr
   T = 150.0_pr

   call model%pressure(n, V, T, P)
end program main
```

A related complete implementation is available in
`example/extra/adiff_pr76.f90`. The example here makes the total-volume
initializer explicit; compare conventions rather than copying a molar
initializer unchanged into a new extensive model.

@note
The `Ar` function must be a function returning `type(hyperdual)`, not a
subroutine with an output argument. The return value is the residual Helmholtz
energy itself, not a derivative array.
@endnote

## Tapenade source-to-source AD

Tapenade reads a conventional Fortran routine and writes new Fortran routines
that evaluate derivatives. Unlike hyperdual numbers, the derivative operations
are explicit in the generated source. This can improve runtime performance,
but requires a separate generation step and review of the generated code.
Tapenade must be installed and available as `tapenade` in `PATH`; see the
[Tapenade download and user documentation](https://tapenade.gitlabpages.inria.fr/userdoc/build/html/download.html).

### What the `yaeos` wrapper needs

The abstract type [[ArModelTapenade]] adapts generated routines to the standard
`ArModel` derivative interface. A concrete model supplies the primal routine
`ar` and four derivative routines:

| Binding | Role in the wrapper |
|:--|:--|
| `ar` | Evaluate the residual Helmholtz energy \(A^r\). |
| `ar_d` | Forward directional derivatives. |
| `ar_b` | Reverse-mode derivatives with respect to the input variables. |
| `ar_d_b` | Mixed forward/reverse derivatives, used for composition Hessians and mixed derivatives. |
| `ar_d_d` | Second forward directional derivatives. |

The routines have strict interfaces; use the library API and the working PC-SAFT
implementation as references rather than guessing argument order:
`src/adiff/autodiff_api/tapenade_ar_api.f90` and
`src/models/residual_helmholtz/saft/pc_saft.f90`.
The derived type must also implement the inherited `get_v0` initializer, which
provides a liquid-volume starting estimate to the volume solver.

### Write and prepare the primal routine

Start from `tools/tapenade_diff/template.f90` in the repository.
The template is a scaffold, not a compilable model: replace its placeholders,
implement the model parameters and constructor, and write the primal routine
with the expected argument names and order:

```fortran
subroutine ar(model, n, v, t, arval)
   type(ModelName), intent(in) :: model
   real(8), intent(in) :: n(:)
   real(8), intent(in) :: v
   real(8), intent(in) :: t
   real(8), intent(out) :: arval
   ! Calculate A^r(n, v, t) here.
end subroutine ar
```

The example uses `real(8)` because the generator workflow expects that kind in
its input. Its post-processing script converts generated declarations to
`real(pr)`. Keep the interface names unchanged: the script and Tapenade
directives refer to them. The residual Helmholtz routine should compute the
primal value only; do not hand-code derivatives in `ar`.

For derived-type parameters, the template recommends copying values to local
variables with explicit shapes based on `size(n)`. This gives Tapenade clear
shape information during differentiation. Keep the routine deterministic and
avoid I/O, hidden state changes, and unsupported external calls in the
differentiated calculation.

### Generate the derivative source

Inspect the script **before running it**, in a disposable working copy. Its
paths are relative to the current directory. The intended invocation from the
generator directory is:

```bash
cd tools/tapenade_diff
bash gen_tapemodel.sh /absolute/path/to/your_model.f90
```

There is a path caveat in the checked-in script: `ADFirstAidKitDIR` is
`../src/adiff/autodiff_api/tapenade`, but from `tools/tapenade_diff/` the source
tree is at `../../src/adiff/autodiff_api/tapenade`. Correct the helper path in
your working copy (or otherwise supply the expected layout) before invoking
the generator. The command above is not a promise of an unattended build.

The script performs a sequence of tangent and reverse Tapenade passes and
post-processes the result into `tapeout/your_model_diff.f90`. **It clears files
in its `tapeout/` directory before generating output.** Save any work there
elsewhere before running it. Check the script and the installed Tapenade version
before use; this is a developer workflow, not a command to run on an existing
project output directory. The script also removes intermediate sources and
Tapenade message files; preserve diagnostics separately when investigating
generation problems. It does not use fail-fast shell handling, so check every
pass for errors rather than treating the final file's existence as success.

The current pass sequence is tangent, nested tangent, third tangent,
reverse-of-tangent, then reverse of the primal. The final file is renamed and
post-processed: module flags are stripped, `REAL*8`/kind suffixes are adjusted,
derived-type dummy declarations are changed toward polymorphic `class`, and
external declarations are removed. These textual substitutions are not a
Fortran interface checker; inspect the result, including `use` statements and
the availability of the Tapenade support interfaces.

### Review and integrate the result

Generation is not the end of the process. Review the output against the wrapper
interfaces, then make the necessary project-specific fixes:

1. Put the generated routines in a module that can access the model type and
   `pr` kind.
2. Bind `ar`, `ar_d`, `ar_b`, `ar_d_b`, and `ar_d_d` to the derived type.
3. Check that `model` is the passed-object argument and that all other dummy
   arguments match the deferred interfaces exactly.
4. Check array extents, kinds, and `intent` declarations; ensure the generated
   code uses the correct component count and does not retain temporary
   Tapenade-specific names or flags.
5. Add the required `get_v0` implementation and compile the model with the
   library.

`tools/tapenade_diff/tapeout/pr_diff.f90` is a post-processed Peng-Robinson
example. The production PC-SAFT implementation
is another reference, but its generated derivatives are integrated into a
larger module and should not be treated as raw Tapenade output.

### Seeds, adjoints, and repeated calls

Normally call `residual_helmholtz` and let the adapter manage these details.
When testing generated callbacks directly, follow
`src/adiff/autodiff_api/tapenade_ar_api.f90` exactly. The generated reference
implementation in `tools/tapenade_diff/tapeout/pr_diff.f90` shows the
conventions concretely.

A tangent seed selects a direction. Zero all tangent inputs, then set, for
example, `vd = 1` to obtain the volume derivative in `arvald`.

Reverse sweeps consume their output seed. In the generated `ar_b`, the
incoming energy seed `arvalb` is set to zero before the routine returns, and
in `ar_d_b` both `arvalb` and `arvaldb` are zeroed. A seed value read after
the call is therefore not the value you passed in.

Do not rely on how a particular generated routine initializes its input
adjoints. In the checked-in `ar_b` routines, `nb` and `tb` are first assigned
and then incremented within the routine (for example `tb = tb + SUM(...)`),
but other generated code or a hand-edited copy may accumulate onto the incoming
value instead. The safe rule is the one the wrapper follows: clear `nb`, `vb`, and `tb` before **each** independent reverse sweep,
and set `arvalb = 1` again before every `ar_b` call.

For `ar_d_b`, clear both adjoint families and all unused seeds. Set a tangent
direction and `arvaldb = 1`, leaving `arvalb = 0`. With `nd(i) = 1`, the
wrapper reads a composition-Hessian row from `nb`; `ndb` supplies the mole
gradient. Volume and temperature adjoints provide the corresponding mixed
derivatives. `arvalb` and `arvaldb` are not interchangeable.

Nested forward calls use two seed directions: `vd = vd0 = 1` gives `ArV2`;
`vd = 1` and `td0 = 1` gives `ArTV`, with all other seeds zero. Read the
second derivative from `arvaldd`, not from a first-tangent output.

Do not assume a reverse sweep returns the primal energy. The wrapper calls
the primal routine separately when `Ar` is requested. All five Tapenade
bindings are **subroutines**, unlike the hyperdual `Ar` function.

The adapter's `reset_vars` clears tangent, adjoint, mixed, and scalar workspace
before independent sweeps. Preserve that behavior in custom wrappers. Nonzero
stale adjoints can make a first call look correct but repeated calls or requests
for several derivatives fail. Test repeated calls and combined requests, not
only a single isolated Hessian evaluation.

### Validate the derivatives

Compare generated derivatives with `numeric_ar_derivatives` at multiple
compositions, temperatures, and volumes; see the analytic guide for its required
energy output and finite-difference step choices. Thermodynamic identities are
a separate check. In a program with a declared scalar `eq31`, for example:

```fortran
use yaeos__consistency_armodel, only: ar_consistency

call ar_consistency(model, n, V, T, eq31=eq31)
print *, "Fugacity/Gibbs identity residual:", eq31
```

Collect all five identity outputs and check independent optional calls as shown
in the analytic guide. Small identity residuals do not by themselves prove that
the derivatives match the intended primal equation or experimental data.

Also test limiting cases (pure-component compositions, dilute components, and
volumes near the model's lower limit) and verify that the model returns finite
pressures and properties before using it in flash or envelope calculations.

Turn these checks into a permanent regression test once they pass. The
[testing page of the developer handbook](../../developer/testing.html) shows how
the existing models are tested, which tolerances the suite uses, and how to
choose a step size that does not hide a wrong derivative.

### Choosing between the two approaches

Hyperdual AD is usually simpler for prototyping and keeps the derivative
implementation close to the primal equations. Tapenade takes more setup and
maintenance, but can generate efficient routines for production models. Both
approaches need validation; neither makes an incorrect primal equation
thermodynamically correct.
