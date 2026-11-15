---
title: Mixing Rules
---

[TOC]

# Turning pure-fluid parameters into a mixture

A cubic EOS first calculates each pure component's attraction and covolume.
A mixing rule answers the next question: how do unlike molecules interact?
It combines those parameters into a single-fluid representation of the mixture.
Changing it can change phase behavior even if every pure-fluid curve stays
unchanged. It is not the same as an [alpha function](alpha.html), which sets the
temperature dependence of each pure component's attraction.

Every [[CubicEoS]] contains a [[CubicMixRule]]. Start with quadratic mixing
for a simple baseline. Consider excess-Gibbs mixing when liquid non-ideality
and an appropriate activity-coefficient parameterization matter. Neither choice
replaces checking against experimental mixture data.

## Which rules are available?

| Rule | Unlike-interaction data | What changes relative to QMR? |
|------|-------------------------|-------------------------------|
| [[QMR]] | Two dimensionless `(nc,nc)` matrices `k`, `l` | The usual pairwise baseline |
| [[QMRTD]] | `k`, `l`, `k0`, `Tref`, all `(nc,nc)` | Attraction corrections vary with temperature |
| [[CMR]] | Two dimensionless `(nc,nc,nc)` tensors `k`, `l` | Composition enters through component triplets |
| [[CMRTD]] | `k`, `l`, `k0`, `Tref`, all `(nc,nc,nc)` | Temperature-dependent triplet corrections |
| [[HV]] | A `GeModel` plus the cubic's `b` and `del1` | Matches an infinite-pressure excess-energy relation |
| [[MHV]] | A `GeModel`, the cubic's `b`, a nonzero `q`, and optional `lij` | Uses an approximate zero-pressure relation |
| [[HV_NRTL]] | Specialized NRTL-HV parameters plus a pair-selection mask and `kij` | Allows selected pairs to use QMR-style corrections |
| [[sDDLC]] | QMR-temperature-law matrices plus positive segment sizes `q(:)` | The attractive mixture parameter also depends on density |

Here `nc` is the number of components. Temperature scales `Tref` are in K;
NRTL-HV energy parameters have their own units, explained below. These are
alternative formulations, not an accuracy ranking. A fitted parameter belongs
to its EOS, alpha function, mixing rule, and temperature convention together.
In particular, a QMR pair matrix is not a CMR triplet tensor.

## Quadratic mixing (QMR)

Let \(n_i\) be the amount of component \(i\) in mol, \(n=\sum_i n_i\) the total
amount, and \(x_i=n_i/n\) its mole fraction. The pure attraction parameter
\(a_i(T)\) has units bar L²/mol² and covolume \(b_i\) has units L/mol.
The dimensionless binary corrections \(k_{ij}\) and \(l_{ij}\) modify attraction
and covolume, respectively. Define the unlike-pair parameters \(a_{ij}\) and
\(b_{ij}\), and the mixture molar parameters \(a_m\) and \(b_m\), by

\[
a_{ij}=\sqrt{a_i(T)a_j(T)}(1-k_{ij}),\qquad
b_{ij}=\frac{b_i+b_j}{2}(1-l_{ij}),
\]
\[
a_m=\sum_i\sum_j x_i x_j a_{ij},\qquad
b_m=\sum_i\sum_j x_i x_j b_{ij}.
\]

Internally the library uses total-volume parameters \(D=n^2a_m\) and \(B=nb_m\):

\[
D=\sum_i\sum_j n_i n_j a_{ij},\qquad
nB=\sum_i\sum_j n_i n_j b_{ij}.
\]

Thus `B` is in L, not L/mol. With zero `lij`, this reduces to
\(B=\sum_i n_i b_i\). Positive `kij` weakens unlike attraction at fixed pure
parameters; it does not mean an equilibrium pressure must change in the same
direction in every calculation.

### Complete example: comparing a binary correction

This methane/n-butane example compares pressure before and after replacing
the default mixing rule. The correction `0.05_pr` is **illustrative**, not a
recommended fitted value. All vectors and both matrix axes use methane first,
n-butane second. Matrices are dimensionless, symmetric and zero on the diagonal.

```fortran
program mixing_comparison
   use yaeos, only: pr, CubicEoS, PengRobinson76, QMR
   implicit none
   type(CubicEoS) :: model
   type(QMR) :: rule
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: kij(2, 2)
   real(pr) :: lij(2, 2)
   real(pr) :: n(2)
   real(pr) :: T
   real(pr) :: V
   real(pr) :: Pbase
   real(pr) :: Pmodified

   tc = [190.564_pr, 425.12_pr]
   pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]
   model = PengRobinson76(tc, pc, w)
   n = [0.3_pr, 0.7_pr] ! mol
   T = 300.0_pr ! K
   V = 10.0_pr ! total L
   call model%pressure(n, V, T, Pbase)

   kij = 0.0_pr
   lij = 0.0_pr
   kij(1, 2) = 0.05_pr
   kij(2, 1) = kij(1, 2)
   rule = QMR(k=kij, l=lij)
   call model%set_mixrule(rule)
   call model%pressure(n, V, T, Pmodified)
   print *, "Default pressure [bar]: ", Pbase
   print *, "Modified pressure [bar]: ", Pmodified
end program mixing_comparison
```

Passing `kij=kij, lij=lij` directly to PR/SRK/RKPR constructors is the simpler
route when no later replacement is needed. For a separately constructed `QMR`,
initialize both matrices; do not leave allocatable fields unset.

As checks, setting the correction to zero should recover the baseline, and
using only one component should recover its pure-fluid prediction. Agreement
at one state does not validate bubble pressures or phase compositions elsewhere.

## Temperature-dependent quadratic mixing (QMRTD)

[[QMRTD]] uses the same covolume rule, but its binary attraction correction
depends on absolute temperature \(T\) (K). Here \(k_{ij}^{\infty}\) is the
high-temperature limiting value, \(k_{ij}^{0}\) is an exponential amplitude
(both dimensionless), and \(T_{ij}^{*}\) is a positive reference temperature in K:

\[
k_{ij}(T)=k_{ij}^{\infty}+k_{ij}^{0}\exp(-T/T_{ij}^{*}).
\]

In particular, the field `k0` is an **amplitude**, not the value at zero
temperature. `k` stores the high-temperature limit. Use the convention in
this equation when transferring fitted parameters from another package.

The following is a **replacement snippet** for the example above: add these
declarations to its specification section, add `QMRTD` to its `use` list,
then replace the construction and assignment of `rule` with these statements.
The numbers are illustrative, not fitted methane/n-butane parameters.

```fortran
! Additional declarations (before executable statements):
type(QMRTD) :: td_rule
real(pr) :: k0(2, 2)
real(pr) :: Tref(2, 2)

! Executable replacement:
k0 = 0.0_pr
k0(1, 2) = 0.02_pr
k0(2, 1) = k0(1, 2)
Tref = 190.0_pr
td_rule = QMRTD(k=kij, l=lij, k0=k0, Tref=Tref)
call model%set_mixrule(td_rule)
```

Keep **all** `Tref` entries positive, including diagonals: the implementation
evaluates an entire matrix expression before selecting pair entries. Symmetry
and component ordering must hold for `k`, `k0`, `l` and `Tref` alike. Avoid
extrapolating temperature fits beyond their data range.

## Cubic composition mixing (CMR and CMRTD)

“Cubic” here describes the **composition mixing polynomial**, not a new
pressure equation. CMR still plugs into `CubicEoS`. Unlike QMR's double sum,
it combines three component indices, including repeated indices:

\[
a_{ijk}=(a_i a_j a_k)^{1/3}(1-k_{ijk}),\qquad
b_{ijk}=\frac{b_i+b_j+b_k}{3}(1-l_{ijk}).
\]
\[
a_m=\sum_{ijk}x_i x_j x_k a_{ijk},\qquad
b_m=\sum_{ijk}x_i x_j x_k b_{ijk}.
\]

The extensive parameters remain `D = n**2 * a_m` and `B = n * b_m`:

\[
D=\frac{\sum_{ijk}n_i n_j n_k a_{ijk}}{n},\qquad
B=\frac{\sum_{ijk}n_i n_j n_k b_{ijk}}{n^2}.
\]

The denominators are essential when amounts do not sum to one. Triplet tensors
are useful even for a **binary** mixture: entries such as `(1,1,2)` and
`(1,2,2)` weight its two unlike composition terms differently. With three
components there are also all-distinct triplets.

Initialize every entry and use the same value for every permutation of a
triplet: for example, `(1,1,2)`, `(1,2,1)`, and `(2,1,1)` must agree. Set the
pure-component entries `(i,i,i)` to zero to preserve the pure-fluid parameters.
With `l=0`, the covolume reduces to linear mixing. With **both** tensors zero,
the CMR attraction is generally **not** QMR: it gives
\((\sum_i x_i a_i^{1/3})^3\), rather than
\((\sum_i x_i\sqrt{a_i})^2\).

CMRTD uses the same exponential convention as QMRTD, but for each triplet:

\[
k_{ijk}(T)=k_{ijk}^{\infty}+k_{ijk}^{0}\exp(-T/T_{ijk}^{*}).
\]

`k0` is again an amplitude, not the correction at zero temperature. Keep
`Tref` positive and respect permutation symmetry in every tensor.

### Complete example: constant and temperature-dependent tensors

The following program constructs both rules. The off-pure entries are
deliberately simple **illustrative parameters**, not a fit to this mixture.

```fortran
program triplet_mixing
   use yaeos, only: pr, CubicEoS, PengRobinson76, CMR, CMRTD
   implicit none
   type(CubicEoS) :: model
   type(CMR) :: constant_rule
   type(CMRTD) :: temperature_rule
   real(pr) :: k3(2, 2, 2)
   real(pr) :: l3(2, 2, 2)
   real(pr) :: amplitude(2, 2, 2)
   real(pr) :: reference_temperature(2, 2, 2)
   real(pr) :: n(2)
   real(pr) :: pressure
   integer :: i
   integer :: j
   integer :: k

   model = PengRobinson76([190.564_pr, 425.12_pr], &
      [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   n = [0.3_pr, 0.7_pr]
   k3 = 0.0_pr
   l3 = 0.0_pr
   amplitude = 0.0_pr
   reference_temperature = 200.0_pr
   do i = 1, 2
      do j = 1, 2
         do k = 1, 2
            if (i /= j .or. j /= k) then
               k3(i, j, k) = 0.02_pr
               amplitude(i, j, k) = 0.01_pr
            end if
         end do
      end do
   end do

   constant_rule = CMR(k=k3, l=l3)
   call model%set_mixrule(constant_rule)
   call model%pressure(n, 10.0_pr, 300.0_pr, pressure)
   print *, "CMR pressure [bar]: ", pressure

   temperature_rule = CMRTD(k=k3, l=l3, k0=amplitude, &
      Tref=reference_temperature)
   call model%set_mixrule(temperature_rule)
   call model%pressure(n, 10.0_pr, 300.0_pr, pressure)
   print *, "CMRTD pressure [bar]: ", pressure
end program triplet_mixing
```

Setting `amplitude=0` must make CMRTD reproduce the **CMR** result with the
same tensors. It is not a test against QMR. Also check pure limits and joint
amount/volume scaling before using fitted tensors in equilibrium calculations.

## Mixing with an excess-Gibbs model

An excess-Gibbs model describes how liquid mixing differs from ideal mixing.
Coupling it to a cubic EOS can introduce liquid-mixture information while the
EOS still supplies volumetric properties. A successful coupling needs a
compatible model, parameter convention and component order; arbitrary NRTL
parameters are not automatically transferable between mixing-rule variants.

For the equations below, \(G^E\) is **total** excess Gibbs energy in bar L,
\(R\) is in bar L/(mol K), and \(T\) is in K. The previously defined \(D\),
\(B\), \(n_i\), \(n\), \(a_i\) and \(b_i\) retain their meanings.

### Huron-Vidal

[[HV]] uses an infinite-pressure limiting relation. For a common cubic shape,
let \(\delta_1\) and \(\delta_2\) be its dimensionless denominator parameters,
and define the positive shape factor \(\Lambda\). The sign convention used by
the implementation is

\[
\Lambda=\frac{\ln[(1+\delta_1)/(1+\delta_2)]}{\delta_1-\delta_2},\qquad
\frac{D}{B}=\sum_i n_i\frac{a_i(T)}{b_i}-\frac{G^E}{\Lambda}.
\]

Do not insert a factor `RT` on only one side: every term here has units bar L.
For PR and SRK, all components share the same cubic shape, so this factor is
composition independent. A formulation with composition-dependent shape, such
as RKPR, must also account for its derivatives; the displayed relation alone
is not a complete derivative implementation.

### Modified Huron-Vidal (MHV)

[[MHV]] uses an approximate zero-pressure relation. Its dimensionless constant
\(q\) belongs to the chosen cubic/mixing formulation, not to a freely fitted
binary interaction. With the library's convention,

\[
\frac{D}{RTB}=\sum_i n_i\frac{a_i(T)}{RTb_i}
+\frac{1}{q}\left[\frac{G^E}{RT}
+\sum_i n_i\ln\left(\frac{B}{nb_i}\right)\right].
\]

The constructor is `MHV(ge=ge_model, b=model%b, q=q)` with optional `lij`.
Build `ge_model` with parameters appropriate to the intended formulation, then
use `call model%set_mixrule(mixrule)`. See the complete repository example
`example/basics/7_huron_vidal.f90` for setup. [[PSRK]] configures MHV with
`q=-0.64663_pr` and PSRK-UNIFAC itself.

`MHV` supports a repulsive `lij` matrix in the implementation, but the
standard MHV derivation uses linear covolume mixing. Keep `lij=0` unless you
have a justified, consistently parameterized extension. `q` must be nonzero;
its appropriate value depends on the cubic and the chosen matching convention.

### Specialized NRTL-HV coupling

[[HV_NRTL]] is a mixing rule; [[NRTLHV]] is the excess-Gibbs model it contains.
They are not interchangeable objects. Its constructor requires all of
`b`, `del1`, `alpha`, `gji0`, `gjiT`, `use_kij`, and `kij`.

For a pair with `use_kij(j,i)=.false.`, the rule uses the supplied NRTL-HV
energy differences and non-randomness parameter. For a selected pair it
constructs the QMR-style interaction internally from the cubic attraction,
covolumes, and `kij`, and sets the corresponding non-randomness parameter to
zero. Supply both matrix directions consistently for each selected pair.
Initialize `kij` even if no pair selects it.

The NRTL-HV law uses referenced differences
\(\Delta g_{ji}=g^0_{ji}+g^T_{ji}T\), with `gji0` in bar L/mol and `gjiT`
in bar L/(mol K). These are **not** the ordinary NRTL constructor's
dimensionless `a` and kelvin-valued `b` matrices. The
[NRTL-HV guide](../../excessmodels/nrtl-hv.html) explains the conventions.

### Complete example: installing each excess-Gibbs rule

This example uses SRK and illustrative liquid-interaction parameters to show
the three setup paths. It does not imply that HV, MHV, and NRTL-HV with these
different parameter choices should predict the same pressure or phase diagram.

```fortran
program excess_gibbs_mixing
   use yaeos, only: pr, CubicEoS, SoaveRedlichKwong, NRTL, HV, MHV, HV_NRTL
   implicit none
   type(CubicEoS) :: model
   type(NRTL) :: ge
   type(HV) :: hv_rule
   type(MHV) :: mhv_rule
   type(HV_NRTL) :: specialized_rule
   real(pr) :: a(2, 2)
   real(pr) :: b_temperature(2, 2)
   real(pr) :: alpha(2, 2)
   real(pr) :: g0(2, 2)
   real(pr) :: gT(2, 2)
   real(pr) :: kij(2, 2)
   real(pr) :: lij(2, 2)
   logical :: use_kij(2, 2)
   real(pr) :: n(2)
   real(pr) :: pressure

   model = SoaveRedlichKwong([190.564_pr, 425.12_pr], &
      [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   n = [0.3_pr, 0.7_pr]
   a = 0.0_pr
   a(1, 2) = 0.2_pr
   a(2, 1) = 0.4_pr
   b_temperature = 0.0_pr
   alpha = 0.3_pr
   alpha(1, 1) = 0.0_pr
   alpha(2, 2) = 0.0_pr
   ge = NRTL(a, b_temperature, alpha)

   ! HV uses the public components of its derived-type constructor.
   hv_rule = HV(ge=ge, bi=model%b, del1=model%del1)
   call model%set_mixrule(hv_rule)
   call model%pressure(n, 10.0_pr, 300.0_pr, pressure)
   print *, "HV pressure [bar]: ", pressure

   lij = 0.0_pr
   mhv_rule = MHV(ge=ge, b=model%b, q=-0.593_pr, lij=lij)
   call model%set_mixrule(mhv_rule)
   call model%pressure(n, 10.0_pr, 300.0_pr, pressure)
   print *, "MHV pressure [bar]: ", pressure

   g0 = 0.0_pr
   g0(1, 2) = 10.0_pr
   g0(2, 1) = 20.0_pr
   gT = 0.0_pr
   use_kij = .false.
   kij = 0.0_pr
   specialized_rule = HV_NRTL(b=model%b, del1=model%del1, &
      alpha=alpha, gji0=g0, gjiT=gT, use_kij=use_kij, kij=kij)
   call model%set_mixrule(specialized_rule)
   call model%pressure(n, 10.0_pr, 300.0_pr, pressure)
   print *, "NRTL-HV pressure [bar]: ", pressure
end program excess_gibbs_mixing
```

The cubic supplies `model%b` in L/mol and `model%del1` in its component order.
The separate name `b_temperature` avoids confusing ordinary NRTL's temperature
matrix with those cubic covolumes. The installed rule owns a copy of the data;
modifying a local rule or activity model afterwards does not update the EOS
until you install the changed rule again.

## Density-dependent local compositions (s-DDLC)

[[sDDLC]] lets the attractive mixture parameter depend on volume as well as
temperature and composition. Its picture is that molecular neighbors need not
have the same proportions as the bulk mixture, and that this local weighting
changes with density. It inherits QMRTD's pairwise temperature law and covolume
rule, then adds a positive dimensionless segment-size vector `q(:)`.

The following equations describe the implemented attractive weighting. Define
the total segment amount \(Q=\sum_i n_i q_i\) and segment fractions
\(s_i=n_iq_i/Q\). With the same pair attraction \(a_{ij}(T)\) as QMRTD,

\[
E_{ij}=\exp\left(\frac{a_{ij}Q}{q_iq_jRTV}\right),\qquad
L=\sum_{ij}s_i s_j E_{ij},\qquad
D=\frac{\sum_{ij}n_i n_j E_{ij}a_{ij}}{L}.
\]

`V` is total volume in L. The segment sizes are model parameters, not the
instantaneous number of moles; changing their common scale also changes the
exponential weighting. Use the convention associated with a fitted data set.
At large volume, `E` tends to one and the attraction approaches the corresponding
QMRTD result. In a pure-component limit the weighting cancels from `D`.

Because `D` now depends on volume, its volume and mixed derivatives matter for
pressure and other properties. Use the **`sDDLC` constructor**, which sets the
`is_D_ddlc` flag used by `CubicEoS`; do not merely copy matrices into a different
rule type. Keep `q` and `Tref` positive, and avoid overflows in the exponential
at very high density or outside the fitted range.

### Complete example: install a density-dependent rule

```fortran
program density_dependent_mixing
   use yaeos, only: pr, CubicEoS, PengRobinson76, sDDLC
   implicit none
   type(CubicEoS) :: model
   type(sDDLC) :: rule
   real(pr) :: kij(2, 2)
   real(pr) :: amplitude(2, 2)
   real(pr) :: reference_temperature(2, 2)
   real(pr) :: lij(2, 2)
   real(pr) :: q(2)
   real(pr) :: n(2)
   real(pr) :: pressure

   model = PengRobinson76([190.564_pr, 425.12_pr], &
      [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   kij = 0.0_pr
   amplitude = 0.0_pr
   reference_temperature = 200.0_pr
   lij = 0.0_pr
   q = [1.0_pr, 3.0_pr] ! Illustrative segment sizes, not fitted data.
   n = [0.3_pr, 0.7_pr]
   rule = sDDLC(k=kij, k0=amplitude, Tref=reference_temperature, l=lij, q=q)
   call model%set_mixrule(rule)
   call model%pressure(n, 1.0_pr, 300.0_pr, pressure)
   print *, "s-DDLC pressure [bar]: ", pressure
end program density_dependent_mixing
```

This is a homogeneous-state property example, not evidence of a stable phase
split. Compare with the base rule at large volume, then check volume and
composition derivatives before moving to a saturation or envelope calculation.

## Extending and checking a rule

To change the attractive combining rule, extend `QMR` and override its
type-bound `aij` procedure. `get_aij` is an interface name, **not** an assignable
procedure-pointer component. The override must supply pair parameters and their
first and second absolute-temperature derivatives; see [[QMR]] and the source
interface. Incorrect derivatives can leave pressure plausible while breaking
heat capacities or equilibrium Jacobians.

For QMR-family corrections, use symmetric matrices with zero self-interactions
even where the implementation only reads one triangle; CMR-family tensors
need permutation symmetry. This does **not** require the directed interaction
matrices of an activity model to be symmetric. Verify the corresponding
pure-fluid and zero-correction limits,
and check pressure and derivative consistency before attempting a difficult
flash. A solver's convergence says nothing about whether fitted interactions
remain physically adequate outside their calibration range.

For a CMR extension, override `aijk` and return rank-three parameter and
temperature-derivative tensors; do not use QMR's rank-two interface. For an
entirely different mixing expression, extend `CubicMixRule` and implement
`Dmix`, `Bmix`, and `D1mix`. Density dependence also requires the density-aware
energy path. The [Fortran developer guide](../../../developer/fortran.html)
describes these contracts.

## Next reading

Read [alpha functions](alpha.html) for the pure-component part,
[excess-Gibbs models](../../excessmodels/index.html) for liquid non-ideality,
[phase equilibrium](../../phase_equilibrium/index.html) for using the model
in saturation and flash calculations, and
[thermodynamic fundamentals](../../fundamentals.html) for the shared
conventions. The cubic theory reference is Michelsen and Mollerup (2007),
listed on the [cubic overview](index.html).
