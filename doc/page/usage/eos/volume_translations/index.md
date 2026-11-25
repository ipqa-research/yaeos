---
title: Volume translations
ordered_subpage: constant.md
                 linear.md
---

[TOC]

# Volume-translated models

An equation of state may describe phase equilibria well but predict liquid
densities less accurately. A volume translation changes the volume associated
with a state of an existing EOS and can be fitted or estimated to improve its
volumetric predictions. The method is not limited to cubic equations; the
transformation below applies to a residual Helmholtz model that is valid at the
corresponding base-model volume [1].

Available translations include [constant](constant.html) and
[linear temperature-dependent](linear.html) component shifts. The Helmholtz
transformation below follows Hammer's SINTEF memo [1].

## Volumes, amounts and sign convention

Let \(\mathbf n=(n_1,\ldots,n_m)\) be component amounts in mol and
\(N=\sum_i n_i\) the total amount. The physical total volume is \(V\), in L.
The base EOS is evaluated at a different total volume \(W\):

\[
W=V+C(\mathbf n,T),\qquad V=W-C(\mathbf n,T).
\]

The extensive displacement \(C\) is in L and must scale linearly when all
amounts are multiplied by the same factor. At fixed \(T\), \(P\), and
composition, a positive \(C\) makes the physical volume smaller than the
corresponding base volume. Both \(V\) and \(W\) must lie in the physical domain
of the transformation and the base EOS.

The translation \(C\) is independent of volume and pressure in this
formulation. It may depend on \(T\) and on the amounts. The constant model
uses \(C=\sum_i n_i c_i\), with fixed molar shifts \(c_i\) in L/mol.

## From pressure to residual Helmholtz energy

At fixed \(T,\mathbf n\), \(dW=dV\), so the intended translated pressure is

\[
P_{\rm tr}(\mathbf n,V,T)=P_{\rm base}(\mathbf n,V+C,T).
\]

yaeos models instead supply the extensive residual Helmholtz energy \(A^r\)
in bar L. Using the ideal-gas reference at the **physical** volume \(V\) gives
[1, section 2]

\[
\boxed{A^r_{\rm tr}(\mathbf n,V,T)=
A^r_{\rm base}(\mathbf n,W,T)+NRT\ln\!\left(\frac{V}{W}\right).}
\]

The logarithm is essential: substituting \(W\) into the base residual energy
alone would leave the ideal reference at the wrong volume. Differentiating
the expression above and using
\(P=NRT/V-\partial A^r/\partial V\) recovers the intended pressure.
When \(C=0\), the translated model reduces to the base model.

[[VolumeTranslatedModel]] implements this transformation and extends
[[ArModel]]. A concrete translation supplies \(C\) and its derivatives; the
inherited methods then provide pressure, fugacity coefficients, volume and
residual caloric properties. At fixed \(P,T,\mathbf n\), its `volume` method
delegates to the base model's solver and subtracts \(C\). It does not perform
an independent search for another translated root.

For the detailed derivatives, see Hammer's memo [1, sections 2 and 4]. It
uses the reduced residual Helmholtz energy \(F=A^r/(RT)\), whereas yaeos
returns the extensive residual energy \(A^r\).

## Equilibrium and limitations

For the **linear** rule \(C=\sum_i n_i c_i\), all phases at a common
\(T,P\) receive the same component-specific fugacity correction
\(\ln\phi_{i,\rm tr}=\ln\phi_{i,\rm base}-Pc_i/(RT)\). It cancels in
isofugacity equations, so corresponding equilibrium compositions and
conditions at fixed \(T,P\) are unchanged when the required roots are
physically admissible [2, 3]. This is not a general claim for nonlinear
composition-dependent translations [3], and translated phase volumes do
change.

The parameterization and its range matter. Constant shifts cannot in general
remove the near-critical density discrepancy [2]. Temperature-dependent
translations can alter caloric properties and may cause unphysical behavior,
including crossing isotherms and negative heat capacities, for some
parameterizations and conditions [4].

## References

1. M. Hammer, *Volume shift for generic EOS*, SINTEF Energy Research memo,
   15 June 2021, sections 2 and 4. See also
   [SINTEF - Thermopack](https://github.com/thermotools/thermopack).
2. A. Péneloux, E. Rauzy, R. Fréze, “A consistent correction for
   Redlich-Kwong-Soave volumes,” *Fluid Phase Equilibria* 8 (1982), 7–23.
   [doi:10.1016/0378-3812(82)80002-2](https://doi.org/10.1016/0378-3812(82)80002-2).
3. R. Privat, J.-N. Jaubert, Y. Le Guennec, “Incorporation of a volume
   translation in an equation of state for fluid mixtures,” *Fluid Phase
   Equilibria* 427 (2016), 414–420.
   [doi:10.1016/j.fluid.2016.07.035](https://doi.org/10.1016/j.fluid.2016.07.035).
4. J. Ahlers, J. Gmehling, “Development of an universal group contribution
   equation of state I,” *Fluid Phase Equilibria* 191 (2001), 177–188,
   especially sections 2–3.
