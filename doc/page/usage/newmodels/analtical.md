---
title: Analytical derivatives
---

[TOC]

# Adding a model with analytical derivatives

The thermodynamic properties in `yaeos` are calculated from the residual
Helmholtz energy \(A^r(\mathbf{n}, V, T)\) and its derivatives. The most
direct way to add a new residual Helmholtz model is to implement those
derivatives analytically. This gives control over shared intermediates and
performance, but requires checking every derivative, including mixed ones.

The other available options are:

- Automatic differentiation with hyperdual numbers, described in
  [Automatic differentiation](autodiff.html).
- Source-to-source differentiation with Tapenade, also described in
  [Automatic differentiation](autodiff.html).

# The `ArModel` contract

A residual Helmholtz model extends the abstract type [[ArModel]] and must
implement two deferred procedures:

| Procedure | Interface | Purpose |
|:----------|:----------|:--------|
| `residual_helmholtz` | [[abs_residual_helmholtz(interface)]] | Calculate \(A^r\) and its derivatives |
| `get_v0` | [[abs_volume_initializer(interface)]] | Provide an initial liquid volume for the volume solver |

```fortran
use yaeos, only: ArModel, pr

type, extends(ArModel) :: MyModel
contains
   procedure :: residual_helmholtz => my_ar
   procedure :: get_v0 => my_v0
end type MyModel
```

The `residual_helmholtz` procedure receives the mole numbers `n`, the total
volume `V`, and the temperature `T`. All the outputs are optional. The following
derivatives can be requested:

| Output | Meaning |
|:-------|:--------|
| `Ar` | \(A^r(\mathbf{n}, V, T)\) |
| `ArV` | \(\left(\partial A^r/\partial V\right)_{T,\mathbf{n}}\) |
| `ArT` | \(\left(\partial A^r/\partial T\right)_{V,\mathbf{n}}\) |
| `ArV2` | \(\left(\partial^2 A^r/\partial V^2\right)_{T,\mathbf{n}}\) |
| `ArT2` | \(\left(\partial^2 A^r/\partial T^2\right)_{V,\mathbf{n}}\) |
| `ArTV` | \(\partial^2 A^r/\partial T \partial V\) |
| `Arn(i)` | \(\left(\partial A^r/\partial n_i\right)_{V,T}\) |
| `ArVn(i)` | \(\left(\partial^2 A^r/\partial V \partial n_i\right)_{T}\) |
| `ArTn(i)` | \(\left(\partial^2 A^r/\partial T \partial n_i\right)_{V}\) |
| `Arn2(i,j)` | \(\left(\partial^2 A^r/\partial n_i \partial n_j\right)_{V,T}\) |

The implementation should calculate only the requested outputs. A common
pattern is:

```fortran
if (present(Ar)) Ar = ...
if (present(ArV)) ArV = ...
! ...
```

The `get_v0` function supplies a volume scale that the library's solvers use
to find the liquid root. Read `src/models/solvers/volume.f90` and
`src/models/residual_helmholtz/ar_models.f90` before choosing a value,
because the two solver families use it differently:

- The generic `ArModel%volume` routine uses `get_v0` as the **initial guess**
  for a Newton solve on volume, multiplied by a small factor. A value near the
  dense-liquid volume works well.
- `volume_michelsen` uses the returned value as a **covolume** `B` in the
  reduced coordinate `zeta = B/V`, requires `0 < B < V`, and uses `B` to bound
  the search interval.

In both cases the value must be an **extensive** volume in liters: it has to
grow with the amount of material, because the solver's volume is a total
volume. The cubic models return the mixture covolume computed by their mixing
rule, and PC-SAFT returns the close-packed segment volume
(`src/models/residual_helmholtz/saft/pc_saft.f90`). A model whose parameters
are defined per mole must convert them to a total quantity before returning.

For a model with a per-component covolume `b`, the total covolume is the
natural choice:

```fortran
function my_v0(self, n, P, T) result(v0)
   class(MyModel), intent(in) :: self
   real(pr), intent(in) :: n(:)
   real(pr), intent(in) :: P
   real(pr), intent(in) :: T
   real(pr) :: v0
   v0 = sum(n * self%b)
end function my_v0
```

Returning zero, a negative value, or a quantity that does not scale with `n`
does not fulfill the contract. A molar value is the classic mistake: it looks
right on a one-mole basis and then misleads the solver for any other amount.
A wrong `get_v0` shows up as a solver that fails to converge or that converges
to the wrong branch, not as a wrong pressure, so test it directly: the volume
found for a doubled amount at the same density should also double, and both a
liquid and a vapor root should be reachable at the same state.

## Units and differentiation variables

`n` contains mole amounts [mol], `V` is total volume [L], and `T` is in K.
`Ar` is in bar L; its derivatives carry those units divided by the units of
the differentiated variables. Thus `ArV` is in bar, `ArT` in bar L/K,
`Arn` in bar L/mol, `ArV2` in bar/L, and `Arn2` in bar L/mol².
The array extents are `size(n)` or `size(n), size(n)` as appropriate.

All mole derivatives hold the other mole amounts fixed. If you write the
equation using mole fractions, differentiate `x = n/sum(n)` too. Likewise,
using molar volume internally requires differentiating `v = V/sum(n)`.
Never normalize the input amounts and silently change the meaning of the
derivatives.

## Implementation sequence

1. Derive the extensive energy from the intended EOS and confirm the ideal-gas
   reference. List all composition and temperature dependencies of parameters.
2. Create the concrete type and constructor. Allocate parameter arrays, populate
   the inherited `components` object, and check input sizes and physical ranges.
3. Implement `Ar` first, then first derivatives, then second/mixed derivatives.
   Use local intermediates rather than optional outputs as scratch storage.
4. Implement the independent `get_v0` function and test its scaling with `n`.
5. Run the checks below before testing a phase-equilibrium calculation.

# Example: van der Waals equation

The van der Waals equation can be written as:

$$
    P = \frac{nRT}{V - B} - \frac{D}{V^2}
$$

with the residual Helmholtz energy:

$$
    A^r = -nRT \ln\left(1 - \frac{B}{V}\right) - \frac{D}{V}
$$

where:

$$
    n = \sum_i n_i, \qquad
    B = \sum_i n_i b_i, \qquad
    D = \sum_i \sum_j n_i n_j a_{ij}
$$

and the pure-component parameters are:

$$
    a_i = \frac{27}{64}\frac{R^2 T_{c,i}^2}{P_{c,i}}, \qquad
    b_i = \frac{R T_{c,i}}{8 P_{c,i}}
$$

The cross parameter is \(a_{ij} = \sqrt{a_i a_j}(1 - k_{ij})\).
This example assumes **symmetric** `kij`, temperature-independent `aij` and `b`,
positive critical constants, and matching array sizes. The constructor is kept
short for teaching; validate those assumptions in production code. In
particular, `Di = 2*matmul(aij,n)` uses symmetry. The parameters `a` and `aij`
have units bar L²/mol², `b` has units L/mol, `B` is in L, and `D` is in bar L².

The following module implements this model analytically:

```fortran
module vdw_model
   use yaeos, only: pr, R, ArModel, Substances
   implicit none

   type, extends(ArModel) :: VdW
      real(pr), allocatable :: a(:)      !! Pure attractive parameters
      real(pr), allocatable :: b(:)      !! Pure covolumes
      real(pr), allocatable :: aij(:, :) !! Cross attractive parameters
   contains
      procedure :: residual_helmholtz => vdw_ar
      procedure :: get_v0 => vdw_v0
   end type VdW

contains

   type(VdW) function setup_vdw(tc, pc, kij) result(model)
      real(pr), intent(in) :: tc(:) !! Critical temperatures [K]
      real(pr), intent(in) :: pc(:) !! Critical pressures [bar]
      real(pr), optional, intent(in) :: kij(:, :) !! Binary interaction parameters

      integer :: i
      integer :: j
      integer :: nc

      nc = size(tc)

      model%components = Substances(tc=tc, pc=pc)
      model%a = 27.0_pr/64.0_pr * R**2 * tc**2 / pc
      model%b = R * tc / (8.0_pr * pc)

      allocate(model%aij(nc, nc))
      do i = 1, nc
         do j = 1, nc
            model%aij(i, j) = sqrt(model%a(i) * model%a(j))
            if (present(kij)) then
               model%aij(i, j) = model%aij(i, j) * (1.0_pr - kij(i, j))
            end if
         end do
      end do
   end function setup_vdw

   subroutine vdw_ar(&
      self, n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, &
      Arn, ArVn, ArTn, Arn2 &
      )
      class(VdW), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: V
      real(pr), intent(in) :: T

      real(pr), intent(out), optional :: Ar
      real(pr), intent(out), optional :: ArV
      real(pr), intent(out), optional :: ArT
      real(pr), intent(out), optional :: ArTV
      real(pr), intent(out), optional :: ArV2
      real(pr), intent(out), optional :: ArT2
      real(pr), intent(out), optional :: Arn(size(n))
      real(pr), intent(out), optional :: ArVn(size(n))
      real(pr), intent(out), optional :: ArTn(size(n))
      real(pr), intent(out), optional :: Arn2(size(n), size(n))

      real(pr) :: nt
      real(pr) :: B
      real(pr) :: D
      real(pr) :: L
      real(pr) :: Di(size(n))

      nt = sum(n)
      B = sum(n * self%b)
      D = dot_product(n, matmul(self%aij, n))
      L = log(1.0_pr - B/V)
      Di = 2.0_pr * matmul(self%aij, n)

      if (present(Ar)) Ar = -nt*R*T*L - D/V
      if (present(ArV)) ArV = -nt*R*T*B/(V*(V - B)) + D/V**2

      if (present(ArT)) ArT = -nt*R*L
      if (present(ArTV)) ArTV = -nt*R*B/(V*(V - B))

      if (present(ArV2)) then
         ArV2 = nt*R*T*B*(2*V - B)/(V**2*(V - B)**2) - 2*D/V**3
      end if

      if (present(ArT2)) ArT2 = 0.0_pr

      if (present(Arn)) then
         Arn = -R*T*L + nt*R*T*self%b/(V - B) - Di/V
      end if

      if (present(ArVn)) then
         ArVn = -R*T*B/(V*(V - B)) &
            - nt*R*T*self%b/(V - B)**2 + Di/V**2
      end if

      if (present(ArTn)) ArTn = -R*L + nt*R*self%b/(V - B)

      if (present(Arn2)) then
         Arn2 = R*T*(spread(self%b, 1, size(n)) &
            + spread(self%b, 2, size(n)))/(V - B) &
            + nt*R*T*spread(self%b, 1, size(n)) &
            * spread(self%b, 2, size(n))/(V - B)**2 &
            - 2.0_pr*self%aij/V
      end if
   end subroutine vdw_ar

   function vdw_v0(self, n, P, T) result(v0)
      class(VdW), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: P
      real(pr), intent(in) :: T
      real(pr) :: v0

      v0 = sum(n * self%b)
   end function vdw_v0

end module vdw_model
```

With the model defined, all the inherited `ArModel` methods are available:

```fortran
program main
   use yaeos, only: pr
   use vdw_model, only: VdW, setup_vdw
   implicit none

   type(VdW) :: model
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: n(2)
   real(pr) :: V
   real(pr) :: T
   real(pr) :: P

   tc = [190.0_pr, 425.0_pr]
   pc = [46.0_pr, 38.0_pr]

   model = setup_vdw(tc, pc)

   n = [0.4_pr, 0.6_pr]
   V = 1.0_pr
   T = 300.0_pr

   call model%pressure(n, V, T, P)
end program main
```

# Validating a new model

Analytical derivatives are easy to get wrong. The consistency module contains
**different** tools for numerical comparisons, independent optional-output
calls, and thermodynamic identities. Add these imports and declarations to the
example program above, before its executable statements:

```fortran
use yaeos__consistency_armodel, only: ar_consistency, &
   numeric_ar_derivatives, individual_ar_calls

! Declarations (after implicit none and the other declarations):
real(pr) :: Ar_num
real(pr) :: ArV_num
real(pr) :: ArV_exact
real(pr) :: eq31
real(pr) :: eq33(2, 2)
real(pr) :: eq34(2)
real(pr) :: eq36
real(pr) :: eq37
logical :: passed
```

Then, after initializing `model, n, V, T`, run:

```fortran
call model%residual_helmholtz(n, V, T, ArV=ArV_exact)
call numeric_ar_derivatives(model, n, V, T, &
   d_n=1.e-5_pr, d_v=1.e-5_pr, d_t=1.e-2_pr, &
   Ar=Ar_num, ArV=ArV_num)
print *, "ArV: implemented, numerical", ArV_exact, ArV_num

call individual_ar_calls(model, n, V, T, passed)
print *, "Independent optional outputs agree:", passed

call ar_consistency(model, n, V, T, &
   eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37)
print *, "Identity residuals:", eq31, maxval(abs(eq33)), &
   maxval(abs(eq34)), eq36, eq37
```

`numeric_ar_derivatives` uses central differences of the primal energy. `Ar`
is a required output, even when only a derivative is of interest. Request the
other derivatives by keyword and compare each in turn. The example step sizes
are starting choices in mol, L, and K, not universal tolerances. Try several
steps: truncation error dominates for large steps and cancellation for very
small ones. All perturbed states must remain in the equation's domain; do not
let a central mole perturbation produce negative amounts.

`individual_ar_calls` checks whether each output agrees with the same output
requested alongside all the others. It does not check correctness against an
independent equation, and its comparisons use a strict absolute machine-epsilon
threshold. If it fails, inspect actual differences and optional-output handling
before deciding whether this is roundoff or a defect.

`ar_consistency` returns identity residuals that should approach zero; it does
**not** print diagnostics or compare against finite differences. Calling it
without `eqXX` outputs produces no visible report. The outputs test the
fugacity/Gibbs relation (`eq31`), composition-derivative symmetry (`eq33`), a
Gibbs-Duhem sum (`eq34`), and pressure (`eq36`) and temperature (`eq37`)
identities. Their units and scales differ, so use scale-aware tolerances rather
than a single raw threshold. Start at positive pressure away from spinodals:
the helper evaluates fugacity expressions involving `log(Z)` and pressure-based
sensitivities, which can be undefined or ill-conditioned elsewhere.

Finally, verify `Ar(lambda*n, lambda*V, T) = lambda*Ar(n,V,T)` and unchanged
pressure for several positive `lambda`. A one-mole test alone cannot reveal
all missing total-amount factors. Run across compositions and densities; only
then test liquid/vapor volume roots and equilibrium calculations.

## Diagnosing failures

| Symptom | First things to check |
|:--|:--|
| Pressure differs by a constant factor | bar versus Pa, L versus m³, and molar versus total energy |
| Pressure is correct but composition derivatives fail | Differentiation of `sum(n)`, mole fractions, and mixing rules |
| Temperature derivatives fail | Missing temperature dependence of parameters, or differentiating a reduced energy without its `RT` factor |
| Hessian is asymmetric | Index order, asymmetric interaction inputs, and missing mixed terms |
| Outputs fail only when requested alone | Reading absent optional arguments or relying on another output being present |
| NaN near high density | `V <= B`, invalid logarithms, or an unsuitable volume initializer |
| Flash fails after local checks pass | Root selection, initial states, model validity, and conditioning near critical points |

Consistency is necessary but not sufficient: a consistently differentiated
wrong equation can pass identities. Retain independent equation/parameter
checks and comparisons to justified reference data.

# Notes and recommendations

- Follow the derivative naming convention from the
  [style guide](../../contributing/styleguide.html): `ArV`, `ArT`, `ArV2`,
  `ArTV`, `Arn`, `ArVn`, `ArTn`, and `Arn2`.
- Use the `R` constant from `yaeos__constants`; it is consistent with the
  library units (bar, L, K, mol).
- Every model must carry a `Substances` object in `components`. This object
  stores the critical constants, acentric factors, and critical volumes of the
  components.
- If the model has an analytical covolume, implementing `get_v0` as the total
  covolume usually improves the liquid-root volume calculation.
- The optional outputs are independent. Avoid computing derivatives that were
  not requested when the calculation is expensive.
- If the model is intended for phase-equilibrium calculations, all the second
  derivatives listed above are required.
