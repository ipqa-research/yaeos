---
title: Stability analysis
---

[TOC]

# Stability analysis

A homogeneous mixture can satisfy the equation of state and still lower its
Gibbs energy by separating into phases of different compositions. Stability
analysis asks whether creating a small amount of a trial phase would help.
It is the appropriate first question when you do not know whether to run a
flash, and can supply a useful initial composition for that flash.

This is not the same as selecting a `"stable"` volume root at a fixed
composition. Nor does a negative result tell you the final amount of each
phase: a [flash](flash.html) solves that material balance.

`yaeos` provides two routines:

- [[tm]]: the reduced tangent-plane distance of a trial phase with a given
  composition.
- [[min_tpd]]: search from several starting compositions and return the lowest
  trial-phase distance found.

# Tangent-plane distance

The tangent plane describes the chemical potentials of the homogeneous feed.
A trial phase below that plane can reduce the total Gibbs energy. Let \(z_i\)
be normalized feed mole fractions, \(P\) pressure in bar, \(T\) temperature in
K, and \(W_i\) positive trial-phase mole amounts. `tm` evaluates the modified,
dimensionless distance:

$$
    tm(W) = 1 + \sum_i W_i \left(\ln W_i + \ln \hat{\phi}_i(W) - d_i - 1\right)
$$

where

$$
    d_i = \ln z_i + \ln \hat{\phi}_i(z)
$$

Here \(\hat\phi_i\) is a fugacity coefficient and \(d_i\) is the feed's
reduced chemical-potential term. When \(\sum_i W_i=1\), this reduces to
the ordinary normalized tangent-plane distance
\(\sum_i W_i[\ln W_i+\ln\hat\phi_i(W)-d_i]\).

A significantly negative value establishes instability of the homogeneous
reference within the model. A nonnegative numerical search result means
**no destabilizing phase was found**, not proof of a global minimum. The
trivial trial `w=z` has zero distance even far from a phase boundary.

Fragment for an initialized model (place declarations before executable code):

```fortran
use yaeos, only: pr, tm

real(pr) :: z(2)
real(pr) :: w(2)
real(pr) :: P
real(pr) :: T
real(pr) :: tpd

z = [0.4_pr, 0.6_pr]
w = [0.8_pr, 0.2_pr]
P = 50.0_pr
T = 300.0_pr

tpd = tm(model, z, w, P, T)

if (tpd < 0.0_pr) then
   print *, "The feed is unstable"
end if
```

`tm` accepts `GeModel` objects as well as `ArModel` objects. For a `GeModel` the
activity coefficients replace the fugacity coefficients, and pressure has no
effect on the model.

The optional `d` argument allows the caller to provide precomputed values of
\(d_i\). This avoids recalculating the feed fugacity coefficients:

Fragment for an `ArModel`; also declare `Vz`, `lnphi_z(size(z))`, and `d(size(z))`:

```fortran
call model%lnphi_pt(z, T=T, P=P, V=Vz, root_type="stable", lnPhi=lnphi_z)
d = log(z) + lnphi_z

tpd = tm(model, z, w, P, T, d=d)
```

# Finding the minimum

[[min_tpd]] minimizes \(tm\) with respect to the trial-phase mole numbers:

```fortran
subroutine min_tpd(model, z, P, T, mintpd, w, all_minima)
   class(BaseModel), target :: model
   real(pr), intent(in) :: z(:)
   real(pr), intent(in) :: P
   real(pr), intent(in) :: T
   real(pr), intent(out) :: w(:)
   real(pr), intent(out) :: mintpd
   real(pr), optional, intent(out) :: all_minima(:, :)
end subroutine min_tpd
```

The routine returns:

| Output | Meaning |
|:-------|:--------|
| `mintpd` | Minimum tangent-plane distance |
| `w` | Trial-phase composition at the minimum |
| `all_minima` | Array of shape `(nc,nc+1)`: each row holds `nc` composition entries followed by its distance |

`min_tpd` starts once per component with that component enriched in trial
mole amounts, then uses successive substitution. It normalizes each resulting
trial and reports the lowest normalized distance. `all_minima` exposes those
trials, including possible duplicates; they are not guaranteed distinct or
fully converged minima. There is no returned convergence flag. An `(nc,nc)`
array is too small: the implementation writes the distance in column `nc+1`.

```fortran
program stability_example
   use yaeos, only: pr, CubicEoS, PengRobinson76, min_tpd
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(CubicEoS) :: model
   real(pr) :: tc(2)
   real(pr) :: pc(2)
   real(pr) :: w(2)
   real(pr) :: z(2)
   real(pr) :: P
   real(pr) :: T
   real(pr) :: mintpd
   real(pr) :: wtrial(2)
   real(pr) :: minima(2, 3)

   tc = [190.564_pr, 425.12_pr]
   pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]

   model = PengRobinson76(tc, pc, w)

   z = [0.4_pr, 0.6_pr]
   P = 50.0_pr
   T = 300.0_pr

   call min_tpd(model, z, P=P, T=T, mintpd=mintpd, w=wtrial, all_minima=minima)
   if (.not. ieee_is_finite(mintpd)) error stop "Nonfinite stability result"
   if (.not. all(ieee_is_finite(wtrial))) error stop "Nonfinite trial phase"

   print *, "Minimum TPD: ", mintpd
   print *, "Most unstable trial phase: ", wtrial

   if (mintpd < -1.0e-6_pr) then
      print *, "The feed is unstable: a flash should be performed"
   else
      print *, "No significant negative distance found in this search"
   end if
end program stability_example
```

# Using stability analysis to initialize a flash

The normalized `w` is a trial phase, not necessarily a final equilibrium
composition. The following fragment uses it as the initial `y` estimate relative
to feed `z` as `x`. It requires strictly positive `z`; reverse the phase roles
(`k0=z/w`) when necessary for the desired heavy/light assignment. Add imports
for `flash` and `EquilibriumState`, declare `result` and integer `iters` first:

```fortran
call min_tpd(model, z, P=P, T=T, mintpd=mintpd, w=w)

if (mintpd < -1.0e-6_pr) then
   result = flash(model, z, T, p_spec=P, k0=w/z, iters=iters)
end if
```

This is especially useful for:

- Liquid-liquid equilibria, where the Wilson correlation may be a poor initial
  guess.
- Systems near a critical point, where the two phases have similar
  compositions.
- Activity-coefficient models, where the initial K-factors are mandatory.

# Multiphase calculations

[[pt_mp_flash]] uses stability analysis internally: it tests the feed, solves
a two-phase flash, and then tests the resulting phase compositions, adding one
more phase if the reference phase is still unstable. It does not iterate until
global stability is proved, and its no-split path is not usable. The process,
and its limitations, are described on the
[multiphase flash](multiphase.html) page.

# Practical notes

- A truly globally stable feed is not metastable. A feed **reported** as stable
  by a finite search may be metastable if a lower-energy split was missed.
  Inspect the trials and compare independent initializations/model data when
  phase behavior is complicated.
- For an `ArModel`, the trial phase is evaluated with the stable volume root.
  For liquid-liquid calculations, a separate check with the liquid root may be
  necessary.
- The distance is dimensionless. The example's `1e-6` threshold is an
  illustrative numerical screen, not an adjustable solver tolerance or a
  universal physical threshold. The routine's own stopping tests can limit
  resolution near zero.
- Near-zero values may mean a boundary, criticality, convergence to the
  trivial feed trial, or numerical noise. Check trial/feed differences and
  neighboring states before interpreting them.
- Feed and trial logarithms require positive fractions. Remove absent
  components consistently from the model rather than passing exact zeros to
  this implementation.
- `min_tpd` is also available for `GeModel` objects. In that case the pressure
  is only used for the calling convention and does not enter the activity
  coefficient model.

# Next steps and reference

Use an unstable trial to initialize a [flash](flash.html), then verify the split's
fugacity and balance residuals. A negative stability test on either returned
phase can signal a missing phase; see the limitations on the
[multiphase page](multiphase.html).

Michelsen and Mollerup, *Thermodynamic Models: Fundamental and Computational
Aspects*, Tie-Line Publications (2004), is the formulation referenced by the
implementation ([review DOI](http://dx.doi.org/10.1016/j.fluid.2005.11.032)).
