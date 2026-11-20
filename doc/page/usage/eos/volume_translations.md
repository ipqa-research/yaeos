---
title: Volume translations
---

[TOC]

# Translating the volume of an existing EOS

[[ConstantVolumeTranslation]] wraps an existing [[ArModel]] without changing
its parameters. Each component has a constant molar shift \(c_i\), in L/mol,
in the same order as the base model's components. The total shift is

\[
C(n)=\sum_i n_i c_i.
\]

Here `n(:)` contains mole amounts, not necessarily normalized fractions.
The physical volume is \(V=V_{\mathrm{base}}-C\): positive coefficients
reduce the volume at fixed pressure and temperature.

## Complete Fortran example

The shift below is illustrative, not a fitted methane parameter.

```fortran
program constant_translation_example
   use yaeos, only: pr, CubicEoS, PengRobinson76, ConstantVolumeTranslation
   implicit none
   type(CubicEoS) :: base
   type(ConstantVolumeTranslation) :: model
   real(pr) :: V_base
   real(pr) :: V

   base = PengRobinson76([190.564_pr], [45.99_pr], [0.01142_pr])
   model = ConstantVolumeTranslation(base, c=[0.01_pr])

   call base%volume([1._pr], 20._pr, 300._pr, V_base, "vapor")
   call model%volume([1._pr], 20._pr, 300._pr, V, "vapor")
   if (V <= 0._pr) error stop "Translated volume solve failed"
   if (abs(V-(V_base-0.01_pr)) > 1.e-12_pr) error stop "Volume shift mismatch"
   print *, "Base and translated total volumes [L]:", V_base, V
end program constant_translation_example
```

For this one-mole example, the second volume is 0.01 L smaller. The constructor
copies the base EOS, component data, and shift vector. The returned model
extends `ArModel`, so it can be passed to existing property and equilibrium
routines. This constructor is currently exposed through the Fortran API.

## Thermodynamic correction

The correction is applied to residual Helmholtz energy, not just the reported
volume. With \(W=V+C\) and \(N=\sum_i n_i\),

\[
A^r_{\mathrm{translated}}(n,V,T)
=A^r_{\mathrm{base}}(n,W,T)+NRT\ln(V/W).
\]

The parent [[VolumeTranslatedModel]] derives the corresponding properties.
For this constant linear law, \(C_{n_i}=c_i\), while every temperature
derivative and the amount Hessian of `C` are zero.

The `c` vector must match the base component count and order. Negative and zero
coefficients are allowed, but a selected translated volume must remain positive.
If the base solver fails, or its selected root translates to a nonpositive or
nonfinite volume, `volume` returns `-1`; it does not search another root.
Choose coefficients appropriate to the fluid and range: this adapter does not
estimate them or guarantee improved density predictions.
