---
title: Constant volume translation
---

[TOC]

# Constant volume translation

[[ConstantVolumeTranslation]] wraps an existing [[ArModel]] with fixed
component shifts \(c_i\) (not temperature dependant) and the linear combining
rule

\[
C(\mathbf n)=\sum_i n_i c_i
\]

Supply `c(:)` in L/mol in the same order and with the same component count as
the base EOS. Amounts `n(:)` are in mol and need not sum to one; `C` is a
total volume in L. Positive coefficients reduce the physical volume, as
explained in the [general derivation](index.html).

## Fortran example

The following shifts are illustrative, **not fitted parameters**. This
complete program compares a PR76 vapor root with its translated counterpart
on a two-mole basis.

```fortran
program constant_translation_example
   use yaeos, only: pr, CubicEoS, PengRobinson76, ConstantVolumeTranslation
   implicit none
   type(CubicEoS) :: base
   type(ConstantVolumeTranslation) :: model
   real(pr) :: n(2)
   real(pr) :: c(2)
   real(pr) :: V_base
   real(pr) :: V

   base = PengRobinson76(&
      [190.564_pr, 305.32_pr], &   ! Critical Temperatures
      [45.99_pr, 48.72_pr], &      ! Crutucal Pressures
      [0.01142_pr, 0.0995_pr] &    ! Acentric factors
      )

   c = [0.010_pr, 0.020_pr] ! L/mol; illustrative values
   n = [0.8_pr, 1.2_pr]      ! mol
   
   model = ConstantVolumeTranslation(base, c)

   call base%volume(n, 1._pr, 300._pr, V_base, "vapor")
   call model%volume(n, 1._pr, 300._pr, V, "vapor")
   
   if (abs(V-(V_base-dot_product(n,c))) > 1.e-10_pr) error stop "Volume shift mismatch"
   
   print *, "Base and translated total volumes [L]:", V_base, V
end program constant_translation_example
```

## Effect on properties

For corresponding base and translated states at the same \(T,P,\mathbf n\),
with \(W=V+C\), the pressure is the same. The fugacity correction is

\[
\ln\phi_{i,\rm tr}=\ln\phi_{i,\rm base}-\frac{P c_i}{RT}.
\]

Consequently, the linear correction cancels between phases at a common \(T,P\)
and preserves corresponding phase-equilibrium conditions and compositions
when their roots remain admissible [1]. It does **not** leave individual
fugacity coefficients or volumes unchanged. It also shifts residual enthalpy
by \(-PC\), while residual \(C_P\) at corresponding \(T,P,\mathbf n\)
states is unchanged for temperature-independent \(c_i\).

## Estimating component shifts

The constructor **does not implement or choose a correlation**: calculate or
fit each \(c_i\) externally and supply it. The following are different
literature choices, not interchangeable universal defaults. Use \(R\) in
bar L/(mol K), \(T_c\) in K and \(P_c\) in bar to obtain L/mol.

### Soave–Redlich–Kwong: PSRK

For SRK, Horstmann et al. [3, equations 8–9] estimate the constant volume
shift from the component's critical properties:

\[
c_i=0.40768\,\frac{R T_{c,i}}{P_{c,i}}
        (0.29441 - \frac{P_{c,i} v_{c,i}}{R T_{c,i}}).
\]

Here \(v_{c,i}\) is the experimental critical molar volume, so
\(P_{c,i}v_{c,i}/(RT_{c,i})\) is the experimental critical compressibility
factor. Horstmann et al. also mention estimating the shift from experimental
liquid-density data [3].

### Peng–Robinson 1976: Ahlers–Gmehling VTPR correlation

For the temperature-independent PR volume correction, Ahlers and Gmehling [2,
equations 4, 12–13] use the **experimental** critical compressibility
\(Z_{c,i}=P_{c,i}v_{c,i}/(RT_{c,i})\) and obtain

\[
c_i=0.252\,\frac{R T_{c,i}}{P_{c,i}}
              (1.5448Z_{c,i}-0.4024).
\]

The authors fitted and generalized the PR correction against pure-component
saturated-liquid density data; their constant correction targets the
lower-temperature range rather than the immediate critical region. This
expression is not the SRK Péneloux correlation above. Supplying its \(c_i\)
to a PR76 model applies that constant volume correction; it does **not** by
itself construct the full group-contribution VTPR model or its
temperature-dependent T-VTPR variant.

## Practical limits

Choose coefficients and operating conditions appropriate for the components
and the base EOS. A constant correction cannot generally eliminate errors
near a pure-component critical point [1, 2]. The physical volume \(V\) and
the base volume \(W\) must remain positive and within the base EOS domain.
Root selection uses the base model's solver; see the
[translated-model interface](index.html).

## References

1. A. Péneloux, E. Rauzy, R. Fréze, “A consistent correction for
   Redlich-Kwong-Soave volumes,” *Fluid Phase Equilibria* 8 (1982), 7–23.
   [doi:10.1016/0378-3812(82)80002-2](https://doi.org/10.1016/0378-3812(82)80002-2).
2. J. Ahlers, J. Gmehling, “Development of an universal group contribution
   equation of state I: Prediction of liquid densities for pure compounds with
   a volume translated Peng–Robinson equation of state,” *Fluid Phase
   Equilibria* 191 (2001), 177–188.
3. S. Horstmann, A. Jabłoniec, J. Krafczyk, K. Fischer, J. Gmehling,
   “PSRK group contribution equation of state: comprehensive revision and
   extension IV, including critical constants and \(\alpha\)-function
   parameters for 1000 components,” *Fluid Phase Equilibria* 227 (2005),
   157–164. [doi:10.1016/j.fluid.2004.11.002](https://doi.org/10.1016/j.fluid.2004.11.002).
