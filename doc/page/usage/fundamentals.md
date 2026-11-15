---
title: Thermodynamics essentials
---

[TOC]

# Thermodynamics essentials for using yaeos

You do not need to derive an equation of state to use `yaeos`. You do need to
know what you are specifying, what the model returns, and what checks make a
result believable. This page introduces those ideas before the more detailed
model and calculation guides.

If you have not run a calculation yet, start with the
[getting-started example](index.html). You can return here whenever a guide
uses an unfamiliar term.

## A component is not a phase

A **component** is one chemical species in the model, such as methane or
ethane. A **phase** is a physically distinct region, such as a liquid or a
vapor. A liquid mixture can contain several components, and a pure component
can coexist as both liquid and vapor. Two different liquids can also coexist.

For example, a methane–ethane feed has two components whether it is entirely
vapor or has separated into a liquid and a vapor. In the separated state, both
phases generally contain both components, but in different proportions.

This distinction matters when sizing arrays: a two-component composition has
two entries, not one entry per phase.

## Amount, composition, and volume: three common sources of mistakes

### Mole amounts versus mole fractions

Let \(n_i\) be the amount of component \(i\), in moles. The vector of amounts
is \(\mathbf n\), and the total amount is a scalar:

\[
n_{\mathrm{tot}}=\sum_i n_i.
\]

The **mole fraction** is the fraction of that total contributed by a component:

\[
z_i=\frac{n_i}{n_{\mathrm{tot}}},\qquad \sum_i z_i=1.
\]

For `n = [3.0_pr, 7.0_pr]`, the total is 10 mol and the composition is
`z = [0.3_pr, 0.7_pr]`. This is the same composition as
`n = [0.3_pr, 0.7_pr]`, but not the same amount of material.

Property routines such as `pressure(n, V, T, ...)` use amounts. Many equilibrium
entry points instead request a normalized feed composition `z`. Follow the
individual routine's convention rather than assuming all arrays are normalized
automatically. A feed must have a positive total amount and nonnegative entries;
logarithmic formulations may additionally require strictly positive entries.

### Total volume versus molar volume

The symbol \(V\) denotes total volume, in litres. The **molar volume** is
\(v=V/n_{\mathrm{tot}}\), in litres per mole. A property call with 10 mol and
1 L describes a much denser state than a call with 1 mol and 1 L.

To keep the same temperature, pressure, and composition when scaling a system,
multiply every mole amount and the total volume by the same factor:

\[
P(\lambda\mathbf n,\lambda V,T)=P(\mathbf n,V,T).
\]

That is a useful consistency check. When you pass a normalized composition as
the amount vector, its total is 1 mol, so the numerical total volume is also
the numerical molar volume. This convenient one-mole basis does not mean the
two quantities are interchangeable for other amounts.

### Component order must never change silently

If the first constructor entry describes methane and the second ethane, every
composition, parameter vector, and row/column of an interaction matrix must
use that same order. Array sizes alone cannot catch a reversed composition or
an incorrectly transposed interaction matrix.

## Units are part of the input

The default gas constant is
\(R=0.08314462618\;\mathrm{bar\,L\,mol^{-1}\,K^{-1}}\). Consequently:

| Quantity | Default unit | Frequent mistake |
|----------|--------------|------------------|
| Temperature | K | Passing a Celsius value without adding 273.15 |
| Pressure | bar | Passing pascals, MPa, or atmospheres as if they were bar |
| Total volume | L | Passing cubic metres or a molar volume for a multi-mole input |
| Amount | mol | Treating mass fractions as mole fractions |
| Extensive energy | bar·L | Reading the numerical value directly as joules |
| Extensive entropy or heat capacity | bar·L/K | Reading the result directly as J/(mol·K) |

One bar·litre is 100 J. For an extensive energy `energy` and an amount vector
`n`, the molar energy in J/mol is `100.0_pr * energy / sum(n)`. Some equilibrium
result fields are already molar quantities; check their documentation before
dividing again. SAFT segment diameters and model-specific interaction parameters
have their own units, described on their respective model pages.

## What an equation of state does—and does not do

An **equation of state** (EOS) models how a fluid's pressure, temperature,
volume, and composition are related. The ideal-gas equation is the simplest:

\[
P=\frac{n_{\mathrm{tot}}RT}{V}.
\]

Real molecules occupy space and interact. An EOS adds a description of those
effects, making it possible to estimate liquid densities and non-ideal gas
behavior. Different model families make different approximations; adding more
parameters does not automatically make a model better for your mixture.

`yaeos` separates two jobs. A **model** evaluates properties for a proposed
state. An **equilibrium solver** searches for the phases and compositions that
can coexist. Calling `pressure` at an arbitrary volume evaluates the EOS; it
does not ask whether the material would actually stay in a single phase.

### Which model family should I read about first?

| Your starting problem | Starting point | What to check before trusting it |
|-----------------------|----------------|----------------------------------|
| Gas/hydrocarbon properties and vapor–liquid calculations with critical-property data | [Cubic EOS](eos/cubics/index.html) | Suitable component data and binary interactions |
| Natural-gas mixtures using components in the supported set | [GERG2008](eos/gerg2008.html) | Supported species and the model's validated state range |
| Fluids with available segment parameters | [PC-SAFT](eos/pc_saft.html) | Parameter units and which physical contributions are implemented |
| Non-ideal liquid mixtures with fitted pair parameters | [NRTL or UNIQUAC](excessmodels/index.html) | Parameter convention, temperature range, and reference state |
| Liquid-mixture estimates from molecular functional groups | [UNIFAC-like models](excessmodels/unifac_models/index.html) | Group decomposition and available interaction pairs |

These are starting points, not guarantees of suitability. A liquid activity
model by itself does not supply gas pressure or density. For a particular
mixture, compare predictions with reliable measurements whenever possible.

## Why the library talks about energy

The **Helmholtz energy**, \(A\), is a thermodynamic potential with natural
variables amount, volume, and temperature. Its derivatives connect several
properties consistently. `yaeos` EOS models use the **residual Helmholtz
energy** \(A^r\): the difference from an ideal gas at the same
\((\mathbf n,V,T)\). For example,

\[
P=\frac{n_{\mathrm{tot}}RT}{V}
  -\left(\frac{\partial A^r}{\partial V}\right)_{T,\mathbf n}.
\]

You can use `model%pressure` without calculating that derivative yourself.
The energy formulation is why different EOS types can share the same property
and equilibrium interfaces.

**Residual** does not mean an error in a numerical calculation. It means a
departure from an ideal-gas reference. A residual enthalpy is not the total
enthalpy: the matching ideal-gas contribution is needed too. The reference held
fixed—volume or pressure—also matters for some departure properties. The
[property reference](eos/index.html) explains those distinctions in detail.

An **excess Gibbs energy**, \(G^E\), uses a different reference: ideal mixing
of components in specified pure-component reference states. It measures
non-ideal mixing rather than departure from an ideal gas. Do not substitute an
excess quantity for a residual one merely because both describe non-ideality.

## Fugacity and activity: how models compare phases

Think of **fugacity** as a pressure-like measure of a component's tendency to
leave a phase. More precisely, it expresses chemical potential on a chosen
reference scale. At equilibrium a component has the same chemical potential,
and therefore the same fugacity, in all phases in which it is present.

For an EOS phase with composition \(x_i\), pressure \(P\), and fugacity
coefficient \(\phi_i\),

\[
f_i=x_i\phi_i P.
\]

The coefficient is dimensionless and approaches one in the ideal-gas limit.
`lnPhi` contains \(\ln\phi_i\), **not** \(\phi_i\); use `exp(lnPhi)` when
you want the coefficient itself. Logarithms are convenient for comparing
fugacities that can span many orders of magnitude.

For liquid activity models, the analogous correction is the **activity
coefficient** \(\gamma_i\). Under a chosen liquid reference convention,
\(a_i=x_i\gamma_i\) is the component's activity. A value of one means ideal
mixing relative to that reference, not necessarily ideal-gas behavior.

## A volume root is not a phase split

At a specified pressure and temperature, an EOS can have several mathematical
volume solutions. A small-volume root is liquid-like; a large-volume root is
vapor-like. Depending on the state there may be only one physically admissible
root. The requested names `"liquid"` and `"vapor"` select root branches; they
do not guarantee that two different phases exist.

The `"stable"` root requests the model's preferred homogeneous root. A mixture
at that composition may still lower its Gibbs energy by separating into phases
with different compositions. That is the question answered by
[stability analysis](phase_equilibrium/stability.html) and
[flash calculations](phase_equilibrium/flash.html), not by a volume call alone.

## Match the calculation to the question

| Question | Calculation | Main result |
|----------|-------------|-------------|
| What is the pressure of this proposed homogeneous state? | Property evaluation | Pressure at the supplied amount, volume, and temperature |
| At what pressure does the first bubble or drop appear? | Saturation point | A phase boundary and an incipient phase composition |
| How much liquid and vapor form at this pressure and temperature? | Flash | Phase compositions and fractions |
| How does that boundary change as conditions change? | Phase envelope | A sequence of boundary points |
| Could this homogeneous mixture lower its energy by separating? | Stability analysis | Candidate phases and a stability criterion |
| Where do two phases become indistinguishable? | Critical calculation | A critical state or line |

The [phase-equilibrium overview](phase_equilibrium/index.html) connects these
calculations into a workflow. In a two-phase split the feed composition
\(z_i\), liquid composition \(x_i\), vapor composition \(y_i\), and vapor
mole fraction \(\beta\) satisfy

\[
z_i=(1-\beta)x_i+\beta y_i.
\]

This is a material balance, not a thermodynamic model, and it is one of the
most useful checks on a calculated split. The **K-factor** convention in the
user-facing explanations is \(K_i=y_i/x_i\).

## A checklist before using a prediction

1. Check units, component ordering, and parameter provenance.
2. Check whether the routine expects amounts or normalized compositions.
3. Check its termination/status information before reading result arrays.
4. Check compositions, phase fractions, and material balances.
5. Check that a root or branch corresponds to the problem you intended to solve.
6. Compare with limiting behavior or independent experimental data.

A converged solver has solved its equations to a numerical tolerance. It has
not established that the parameters are appropriate, that a global search found
every possible phase, or that the model is accurate outside its fitted range.
