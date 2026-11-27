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
W=V+C(\mathbf n,T),\qquad V=W-C(\mathbf n,T)
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
P_{\rm tr}(\mathbf n,V,T)=P_{\rm base}(\mathbf n,V+C,T)
\]

yaeos models instead supply the extensive residual Helmholtz energy \(A^r\)
in bar L. To obtain it from the translated pressure, use
\(P=NRT/V-(\partial A^r/\partial V)_{T,\mathbf n}\) and choose the
usual residual reference \(A^r\to 0\) as the volume tends to infinity.
Following Hammer [1, section 2], integration at fixed \(T,\mathbf n\) gives

\[
A^r_{\rm tr}(\mathbf n,V,T)
=\int_V^\infty
 \left[P_{\rm base}(\mathbf n,v+C,T)-\frac{NRT}{v}\right]\,dv
\]

Here \(v\) is a dummy physical volume. Since \(C(\mathbf n,T)\) is
independent of volume, it is constant along this integration path. Add and
subtract the ideal-gas pressure at the corresponding base volume \(v+C\):

\[
\begin{aligned}
A^r_{\rm tr}
={}&\int_V^\infty
 \left[P_{\rm base}(\mathbf n,v+C,T)-\frac{NRT}{v+C}\right]\,dv\\
&+NRT\int_V^\infty\left(\frac{1}{v+C}-\frac{1}{v}\right)\,dv
\end{aligned}
\]

For the first integral, the change of variable \(w=v+C\) gives \(dw=dv\)
and changes the lower limit to \(W=V+C\). Therefore it is the base residual
energy \(A^r_{\rm base}(\mathbf n,W,T)\), with its own ideal-gas reference
at \(W\). The second integral corrects that reference to the **physical**
volume \(V\):

\[
NRT\int_V^\infty\left(\frac{1}{v+C}-\frac{1}{v}\right)\,dv
=NRT\ln\!\left(\frac{V}{W}\right)
\]

Thus, for positive physical and base volumes and a valid integration path,

\[
\boxed{A^r_{\rm tr}(\mathbf n,V,T)=
A^r_{\rm base}(\mathbf n,W,T)+NRT\ln\!\left(\frac{V}{W}\right)}
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

## Helmholtz energy and derivatives

The equations below give the **extensive residual energy** and all its
derivatives implemented by `VolumeTranslatedModel`. They follow the
transformation described in Hammer's memo [1, sections 2 and 4], but use
\(A^r\) in bar L rather than the memo's reduced energy \(A^r/(RT)\).

Throughout, \(W=V+C(\mathbf n,T)\), \(N=\sum_i n_i\), and every
\(A^r_{\rm base}\) and its partial derivatives are evaluated at
\((\mathbf n,W,T)\). A derivative of the base model with respect to \(T\)
or \(n_i\) holds its volume argument \(W\) fixed; derivatives of the
translated model hold its physical volume \(V\) fixed. All other independent
variables are held fixed in each partial derivative. The amounts \(n_i\) are
independent mole amounts, not normalized mole fractions.

### Energy and volume derivatives

\[
A^r_{\rm tr}=A^r_{\rm base}+NRT\ln\!\left(\frac{V}{W}\right)
\]

\[
\frac{\partial A^r_{\rm tr}}{\partial V}
=\frac{\partial A^r_{\rm base}}{\partial W}
+NRT\left(\frac{1}{V}-\frac{1}{W}\right)
\]

\[
\frac{\partial^2 A^r_{\rm tr}}{\partial V^2}
=\frac{\partial^2 A^r_{\rm base}}{\partial W^2}
-NRT\left(\frac{1}{V}-\frac{1}{W}\right)
       \left(\frac{1}{V}+\frac{1}{W}\right)
\]

### Temperature derivatives

\[
\frac{\partial A^r_{\rm tr}}{\partial T}
=\frac{\partial A^r_{\rm base}}{\partial T}
+\left(\frac{\partial A^r_{\rm base}}{\partial W}-\frac{NRT}{W}\right)
 \frac{\partial C}{\partial T}
+NR\ln\!\left(\frac{V}{W}\right)
\]

\[
\begin{aligned}
\frac{\partial^2 A^r_{\rm tr}}{\partial T\,\partial V}
={}&\frac{\partial^2 A^r_{\rm base}}{\partial T\,\partial W}
+NR\left(\frac{1}{V}-\frac{1}{W}\right)\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W^2}
       +\frac{NRT}{W^2}\right)\frac{\partial C}{\partial T}
\end{aligned}
\]

\[
\begin{aligned}
\frac{\partial^2 A^r_{\rm tr}}{\partial T^2}
={}&\frac{\partial^2 A^r_{\rm base}}{\partial T^2}
+2\left(\frac{\partial^2 A^r_{\rm base}}{\partial T\,\partial W}
        -\frac{NR}{W}\right)\frac{\partial C}{\partial T}\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W^2}
        +\frac{NRT}{W^2}\right)
  \left(\frac{\partial C}{\partial T}\right)^2\\
&+\left(\frac{\partial A^r_{\rm base}}{\partial W}
        -\frac{NRT}{W}\right)\frac{\partial^2 C}{\partial T^2}
\end{aligned}
\]

### Mole-amount derivatives

\[
\frac{\partial A^r_{\rm tr}}{\partial n_i}
=\frac{\partial A^r_{\rm base}}{\partial n_i}
+\left(\frac{\partial A^r_{\rm base}}{\partial W}-\frac{NRT}{W}\right)
 \frac{\partial C}{\partial n_i}
+RT\ln\!\left(\frac{V}{W}\right)
\]

\[
\begin{aligned}
\frac{\partial^2 A^r_{\rm tr}}{\partial V\,\partial n_i}
={}&\frac{\partial^2 A^r_{\rm base}}{\partial W\,\partial n_i}
+RT\left(\frac{1}{V}-\frac{1}{W}\right)\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W^2}
        +\frac{NRT}{W^2}\right)\frac{\partial C}{\partial n_i}
\end{aligned}
\]

\[
\begin{aligned}
\frac{\partial^2 A^r_{\rm tr}}{\partial T\,\partial n_i}
={}&\frac{\partial^2 A^r_{\rm base}}{\partial T\,\partial n_i}
+\left(\frac{\partial^2 A^r_{\rm base}}{\partial T\,\partial W}
        -\frac{NR}{W}\right)\frac{\partial C}{\partial n_i}\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W\,\partial n_i}
        -\frac{RT}{W}\right)\frac{\partial C}{\partial T}\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W^2}
        +\frac{NRT}{W^2}\right)
  \frac{\partial C}{\partial T}\frac{\partial C}{\partial n_i}\\
&+\left(\frac{\partial A^r_{\rm base}}{\partial W}
        -\frac{NRT}{W}\right)
  \frac{\partial^2 C}{\partial T\,\partial n_i}
+R\ln\!\left(\frac{V}{W}\right)
\end{aligned}
\]

\[
\begin{aligned}
\frac{\partial^2 A^r_{\rm tr}}{\partial n_i\,\partial n_j}
={}&\frac{\partial^2 A^r_{\rm base}}{\partial n_i\,\partial n_j}\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W\,\partial n_i}
        -\frac{RT}{W}\right)\frac{\partial C}{\partial n_j}\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W\,\partial n_j}
        -\frac{RT}{W}\right)\frac{\partial C}{\partial n_i}\\
&+\left(\frac{\partial^2 A^r_{\rm base}}{\partial W^2}
        +\frac{NRT}{W^2}\right)
  \frac{\partial C}{\partial n_i}\frac{\partial C}{\partial n_j}\\
&+\left(\frac{\partial A^r_{\rm base}}{\partial W}
        -\frac{NRT}{W}\right)
  \frac{\partial^2 C}{\partial n_i\,\partial n_j}
\end{aligned}
\]

The derivatives of \(C\) are supplied by the concrete translation model.
The volume derivatives contain no \(\partial C/\partial V\) terms because
this interface requires \(C\) to be independent of volume.

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
