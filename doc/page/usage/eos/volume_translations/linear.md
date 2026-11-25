---
title: Linear temperature-dependent volume translation
---

[TOC]

# Linear temperature-dependent volume translation

[[LinearVolumeTranslation]] copies a base EOS and uses

\[
c_i(T)=a_i+b_iT,\qquad C(\mathbf n,T)=\sum_i n_i c_i(T).
\]

Provide one intercept \(a_i\) in L/mol and one slope \(b_i\) in L/(mol K)
per component, in base-model order. The temperature is in K. Amounts are in
mol and need not sum to one; \(C\) is a total volume in L. The base EOS is
evaluated at \(W=V+C\), so positive \(C\) reduces the physical volume at
fixed pressure and temperature. Zero slopes recover the
[constant translation](constant.html).

## Example

The parameters below are illustrative, not fitted:

```fortran
program linear_translation_example
   use yaeos, only: pr, CubicEoS, PengRobinson76, LinearVolumeTranslation
   implicit none
   type(CubicEoS) :: base
   type(LinearVolumeTranslation) :: model
   real(pr) :: n(2)
   real(pr) :: a(2)
   real(pr) :: b(2)
   real(pr) :: W
   real(pr) :: V
   real(pr) :: T

   base = PengRobinson76(& 
      [190.564_pr, 305.32_pr], &
      [45.99_pr, 48.72_pr], &
      [0.01142_pr, 0.0995_pr] &
      )
   
   a = [0.01_pr, 0.02_pr] ! L/mol
   b = [-1.e-5_pr, -2.e-5_pr] ! L/(mol K)
   n = [0.8_pr, 1.2_pr] ! mol
   T = 300._pr ! K
   
   model = LinearVolumeTranslation(base, a, b)
   
   call base%volume(n, 1._pr, T, W, "vapor")
   
   call model%volume(n, 1._pr, T, V, "vapor")
   
   if (abs(V-(W-dot_product(n,a+b*T))) > 1.e-10_pr) error stop "Shift mismatch"
end program linear_translation_example
```

## Isotherm crossover and physical consistency

Shi and Li [1, equations 12–16 and section 4.1] study the crossing of
pressure–volume isotherms for **pure-fluid cubic EOS** translations, using
the same sign convention \(v_{\rm tr}=v_{\rm base}-c(T)\). At a given
single-phase state,

\[
\left(\frac{\partial v_{\rm tr}}{\partial T}\right)_P
=\left(\frac{\partial v_{\rm base}}{\partial T}\right)_P-b.
\]

Their no-crossover criterion requires this derivative to be positive over
the temperature and pressure range of interest. If the base EOS has a
positive derivative throughout that range, \(b\leq 0\) is sufficient:
\(b<0\) is the linear case explicitly discussed in section 4.1, while
\(b=0\) is the constant-shift case. A positive \(b\) is not automatically
invalid; the inequality must be checked against the base EOS across the
relevant range. Positive physical and base volumes are also required.

This criterion concerns P–V isotherm crossover, not a guarantee of all
physical properties or of behavior outside the checked range. The paper's
pure-fluid result does not by itself prove the corresponding condition for
mixtures or arbitrary base models. yaeos accepts either sign of `b`; choosing
and checking parameters is the caller's responsibility.

## Reference

1. J. Shi and H. A. Li, “Criterion for determining crossover phenomenon in
   volume-translated equation of states,” *Fluid Phase Equilibria* (2016),
   [doi:10.1016/j.fluid.2016.09.017](https://doi.org/10.1016/j.fluid.2016.09.017),
   accepted manuscript, equations 12–16 and section 4.1.
