---
title: Multiphase flash
---

[TOC]

# Multiphase flash

Some mixtures split into more than two equilibrium phases. Common examples are
water-hydrocarbon systems, gas-condensate systems with a liquid-liquid split,
and systems near a three-phase line. The routine [[pt_mp_flash]] calculates
the equilibrium phase distribution at fixed composition, pressure, and
temperature using stability-based initialization. The current automatic wrapper
has important limitations, described below; it is not an unrestricted global
multiphase-equilibrium solver.

A ternary mixture has three **components**, but need not have three phases.
Conversely, two liquids can each contain all three components. Phase fractions
are amounts of phases on a mole basis, not the fractions of individual
components. A water/hydrocarbon calculation also needs a model and interaction
parameters appropriate for that mixture; a numerical split alone does not
validate cubic-EoS predictions for associating liquids.

# Theory

A multiphase equilibrium satisfies the same isofugacity condition as a
two-phase flash, but for \(N_p\) phases:

$$
   x_{p,i}\,\hat{\phi}_{p,i} = x_{q,i}\,\hat{\phi}_{q,i}
   \qquad (p,q = 1,\dots,N_p,\; i = 1,\dots,N_c)
$$

together with the overall mass balance

$$
   z_i = \sum_{p=1}^{N_p} \beta_p\,x_{p,i}, \qquad
   \sum_{p=1}^{N_p} \beta_p = 1
$$

Here \(N_p\) is the **total** number of phases, \(N_c\) the number of
components, \(x_{p,i}\) the mole fraction of component \(i\) in phase \(p\),
and \(\beta_p\) the fraction of total feed moles in that phase. Each
composition sums to one. \(\hat\phi\) is a fugacity coefficient; equality
of `composition*fugacity coefficient` at common P and T gives equal component
fugacities. These equations define a candidate split, but do not establish
that no other split has lower Gibbs energy.

The implementation stores one phase separately as the **reference phase**
`w`. This is a mathematical reference, not necessarily a zero-amount phase:
its fraction is generally finite after the flash. Internal ratios are
`K(l,i)=x_l(l,i)/w(i)`, not necessarily the two-phase convention `vapor/liquid`.

# `pt_mp_flash`

```fortran
type(MPEquilibriumState) function pt_mp_flash(model, z, P, T)
   class(ArModel), intent(in) :: model
   real(pr), intent(in) :: z(:)   !! Global composition
   real(pr), intent(in) :: P      !! Pressure [bar]
   real(pr), intent(in) :: T      !! Temperature [K]
end function pt_mp_flash
```

## Current wrapper behavior and limitations

The executable code currently:

1. Calls [[min_tpd]] on the feed. It attempts a two-phase solve only if
   `mintpd < -0.001`.
2. Uses one liquid-root main phase and one vapor-root reference phase, solving
   for log composition ratios and log phase fractions.
3. Tests both resulting phase compositions. If the reference-phase test is
   negative beyond the same threshold, it attempts one three-phase solve with
   two liquid-root main phases and a vapor-root reference phase.

It does **not** loop until all phases are stable, despite the workspace's
four-row internal storage. It attempts at most three total phases. The
three-phase initialization also reuses the trial overwritten by the last
stability call, so it requires independent verification.

**Do not use the no-split path as a valid single-phase result.** When the first
stability test is not sufficiently negative, the local `w` and `betas` arrays
are not initialized before the code uses them to assemble the result. That
path therefore has undefined behavior; inspecting `np` afterwards cannot make
it safe. Gate the call with your own
[stability analysis](stability.html) and treat `np < 1` as an error. There is no
public `kind="failed"` or convergence status on `MPEquilibriumState`. The inner
Newton solve has an iteration limit, but its final residual and phase-removal
flags are not exposed by this wrapper.

## What the labels mean

`kinds_x` and `kind_w` record the **volume-root branch** requested from the
equation of state, not a verified physical phase identity. In the example
below, the reference phase is butane-rich even though its requested branch is
labelled `"vapor"`. Use the labels to reproduce the calculation, rather than
inferring phase identity from the label alone. Evaluate the molar volume and
consider the mixture's phase diagram and stability. A pure component's critical
temperature alone does not classify a phase of a mixture.

# Example: gated exploratory call

This complete program checks feed instability first to avoid the known
no-split path. It may deliberately skip the flash. It does not certify that a
returned split is globally stable; use the diagnostics and checks below.

```fortran
program multiphase_example
   use yaeos, only: pr, CubicEoS, PengRobinson76, MPEquilibriumState, pt_mp_flash
   use yaeos, only: min_tpd
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(CubicEoS) :: model
   type(MPEquilibriumState) :: result
   real(pr) :: tc(3)
   real(pr) :: pc(3)
   real(pr) :: w(3)
   real(pr) :: z(3)
   real(pr) :: P
   real(pr) :: T
   real(pr) :: mintpd
   real(pr) :: trial(3)
   real(pr) :: reconstructed(3)
   real(pr) :: lnphi_ref(3)
   real(pr) :: lnphi_phase(3)
   integer :: i

   ! Water, methane, and n-butane
   tc = [647.14_pr, 190.564_pr, 425.12_pr]
   pc = [220.64_pr, 45.99_pr, 37.96_pr]
   w = [0.3443_pr, 0.0115478_pr, 0.200164_pr]

   model = PengRobinson76(tc, pc, w)

   z = [0.1_pr, 0.6_pr, 0.3_pr]
   P = 50.0_pr
   T = 300.0_pr

   call min_tpd(model, z, P, T, mintpd, trial)
   if (.not. ieee_is_finite(mintpd)) error stop "Invalid stability result"
   if (mintpd >= -0.001_pr) then
      print *, "No sufficiently negative TPD; skipping this wrapper."
      stop
   end if

   result = pt_mp_flash(model, z, P, T)
   if (result%np < 1) error stop "No usable split"
   if (.not. allocated(result%betas)) error stop "Missing phase fractions"
   if (size(result%betas) /= result%np+1) error stop "Unexpected fraction count"
   if (.not. all(ieee_is_finite(result%betas))) error stop "Invalid fractions"
   if (.not. all(ieee_is_finite(result%w))) error stop "Invalid reference phase"
   if (.not. all(ieee_is_finite(result%x_l))) error stop "Invalid main phases"
   print *, "Total phase count: ", result%np+1
   print *, "Reference phase: ", result%kind_w, result%w
   print *, "Reference fraction: ", result%betas(result%np+1)

   reconstructed = result%betas(result%np+1)*result%w
   call model%lnphi_pt(result%w, P, T, root_type=result%kind_w, lnPhi=lnphi_ref)
   do i = 1, result%np
      print *, "Phase ", i, " kind: ", result%kinds_x(i)
      print *, "Phase ", i, " composition: ", result%x_l(i, :)
      print *, "Phase ", i, " fraction: ", result%betas(i)
      reconstructed = reconstructed + result%betas(i)*result%x_l(i, :)
      call model%lnphi_pt(result%x_l(i, :), P, T, root_type=result%kinds_x(i), &
         lnPhi=lnphi_phase)
      print *, "Phase ", i, " max log-fugacity mismatch: ", &
         maxval(abs(log(result%x_l(i, :))+lnphi_phase-log(result%w)-lnphi_ref))
   end do
   print *, "Fraction sum: ", sum(result%betas)
   print *, "Maximum component-balance error: ", maxval(abs(reconstructed-z))
end program multiphase_example
```

With the parameters above this run reports `np=2`, that is three phases in
total: two liquid-root main phases and one `"vapor"`-root reference phase.
The phase fractions sum to one and the reconstructed feed matches the input
to about \(10^{-16}\); the log-fugacity mismatches against the reference phase
are of order \(10^{-13}\). Treat those numbers as evidence that the solver
converged to *a* consistent split, not that it found the globally stable one.

# The `MPEquilibriumState` type

The result is a [[MPEquilibriumState]] with:

| Field | Meaning |
|:------|:--------|
| `z` | Global composition |
| `P`, `T` | Pressure [bar] and temperature [K] |
| `np` | Number of main phases, excluding the reference; a valid split has `np+1` total phases |
| `x_l` | Mole fractions of the main phases, `x_l(phase, component)` |
| `w` | Mole fractions of the reference phase |
| `betas` | Mole fractions of the main phases, followed by the reference-phase fraction at index `np+1` |
| `kinds_x` | Kind of each main phase (`"liquid"` or `"vapor"`) |
| `kind_w` | Kind of the reference phase |

The reference phase starts from a stability trial but is a full phase in the
flash result. The main phases have compositions stored in `x_l`. The phase
fractions in `betas(1:np)` refer to the main phases, and `betas(np+1)` is the
reference-phase fraction. The sum of all `np+1` fractions is one.

# Differences from two-phase flash

- [[flash]] is restricted to two phases. A `kind="failed"` flash is not proof
  that the feed is single-phase; it can also mean poor initialization or a
  numerical failure.
- `np=1` means two total phases and `np=2` means three. `np=0` marks the
  no-split branch, which is not a usable single-phase result in the current
  implementation.
- `pt_mp_flash` only supports `ArModel` models. Activity-coefficient models are
  not supported.
- The multiphase flash is more expensive because it performs repeated stability
  analyses and solves a larger system of equations.

# Practical notes

- Check array sizes, finite nonnegative compositions and fractions, composition
  sums, `sum(betas)=1`, and the reconstructed feed. A small balance error alone
  does not show fugacity equality: the reconstruction formula builds balance
  into the solver variables.
- For a two-phase system, `pt_mp_flash` should agree with `flash` within the
  convergence tolerance. Use this as a consistency check.
- Evaluate log fugacities for every phase against the reference using the
  returned root labels. Check that each component's mismatch is small, and
  test the returned phases for additional destabilizing trials.
- Near a three-phase boundary, fractions can approach zero and small changes
  in P or T can change phase count. Scan with small steps and compare against
  independent calculations. This wrapper accepts no previous-state or `k0`
  argument, so it cannot directly reuse a previous split as initialization.
- Stability uses stable volume roots for its trials, but the subsequent flash
  uses the fixed liquid/liquid/vapor root choices described above. It is not
  a general liquid-liquid-only or four-phase automatic solver.
- For phase diagrams with three or more phases, combine `pt_mp_flash` with the
  critical-line and envelope routines to locate the boundaries between regions.
