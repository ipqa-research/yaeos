---
title: UNIFAC Models
ordered_subpage: unifaclv.md
                 dortmund.md
                 psrk.md
---

[TOC]

# UNIFAC Models

## Introduction

UNIFAC predicts liquid nonideality from a molecule's **functional groups** rather
than a separate fitted interaction for every pair of whole molecules. This is
useful for screening solvents or estimating activity coefficients when binary
mixture data are scarce. You supply a chemical decomposition; the constructor
combines it with a built-in table of group sizes, areas, and interactions.

It is not a molecular simulator and does not use a chemical name or atom list
directly. Molecules with the same subgroup counts are indistinguishable to a
given UNIFAC variant, even if their detailed structures differ. A validated
[NRTL](../nrtl.html) or [UNIQUAC](../uniquac.html) fit may be preferable for a
well-characterized binary system; UNIFAC trades that specificity for transfer
of group parameters between compounds.

## Choosing a variant

| Variant | Constructor returning `type(UNIFAC)` | Important difference |
|:--|:--|:--|
| [Original liquid-vapor](unifaclv.html) | `setup_unifac(molecules)` | Original combinatorial term and single-coefficient group interaction law |
| [Dortmund modified](dortmund.html) | `setup_dortmund(molecules)` | Modified size term and temperature-dependent interaction fit |
| [PSRK UNIFAC](psrk.html) | `setup_psrk(molecules)` | Group model/parameter set developed for PSRK EOS use; not itself a full EOS |

The shared Fortran type does **not** make these models interchangeable. Select
the matching group table, interaction set, and temperature function together.
Changing only `d` or copying group IDs from another variant is not a valid model
conversion. Published tables below are useful references; the parameter data
bundled with your installed version determine what is actually available.

Each UNIFAC model is differentiated by the list of available functional groups,
the interaction parameters between functional groups, the value of the exponent
\(d\) in the combinatorial term, and the function of temperature used in the
residual term. Check the section of each model for more specific details of
each model. Let's begin:

[[UNIFAC]] (UNIQUAC Functional-group Activity Coefficients) is an Excess Gibbs
free energy model used to estimate activity coefficients in non-ideal mixtures.
It is particularly useful for predicting the phase behavior of chemical
mixtures, including liquid-liquid equilibrium (LLE) and vapor-liquid
equilibrium (VLE). In this model the Excess Gibbs free energy is calculated
from the contribution of a combinatorial term and a residual term.

$$ \frac{G^E}{RT} = \frac{G^{E,r}}{RT} + \frac{G^{E,c}}{RT} $$

Being:

- Combinatorial: Accounts for the size and shape of the molecules.

- Residual: Accounts for the energy interactions between different functional groups.

"Residual" here names a contribution to **excess** Gibbs energy relative to an
ideal liquid, not a residual departure from an ideal gas as in an EOS.
Differentiating the sum gives \(\ln\gamma_i\). The dimensionless activity is
\(x_i\gamma_i\), and converting to liquid fugacity requires a pure-liquid
reference fugacity. See [Gibbs Excess Models](../index.html) for units, reference
states, and why a liquid activity model alone is not a vapor-liquid flash model.

Each substance of a mixture modeled with [[UNIFAC]] must be represented by a
list a functional groups and other list with the ocurrence of each functional
group on the substance.

### The physical picture

Two ideas explain almost everything about UNIFAC results.

**Size and shape.** A big molecule dilutes itself more than a small one, and a
compact molecule presents a different surface than a floppy one. The
combinatorial term handles this with the two numbers \(r_i\) (relative volume)
and \(q_i\) (relative surface area), which are *sums of subgroup contributions*.
Ethanol and 1-propanol differ in a single CH2 group, so their \(r\) and \(q\)
differ by that group's contributions, and their residual behavior differs only
through the groups they contain.

**Energy of unlike contacts.** Whether an OH group finds water or ethanol
molecules around it changes its chemical potential. The residual term handles
this with group-group interactions, and it is where the "residual" name comes
from. Because interactions are stored between groups rather than molecules,
information transfers: knowing that water likes ethanol's OH group also tells
you something about water with any other alcohol.

The practical consequence is that UNIFAC can rank or screen solvents without
any mixture data, but its accuracy is inherited from the fitted group
interactions. A group pair that was never fitted contributes nothing useful,
and no error message will tell you that.

## Preparing inputs safely

1. Choose the variant before identifying fragments. A **subgroup** identifies a
   specific fragment (for example CH3); a **main group** collects subgroups that
   share interaction parameters. Pass subgroup IDs, not main-group IDs.
2. For each component, set `molecules(i)%groups_ids` and
   `molecules(i)%number_of_groups` to integer arrays of equal length. Counts are
   positive occurrences **per molecule**, not amounts in the mixture. List each
   subgroup once with its total count and do not double-count atoms in compound
   fragments (for example methanol has a dedicated CH3OH group).
3. Check coverage of all required directional main-group interactions, not only
   whether each subgroup exists. A missing interaction is not evidence of an
   ideal or zero-energy contact. The constructors' dimensional checks cannot
   establish scientific coverage or suitability for your mixture.
4. Construct the model once, then vary `n(:)` (mol) and `T` (K). Component order
   in `n` follows `molecules`. Require nonnegative amounts and a positive total.
   Do not substitute mass fractions or Celsius temperatures.

### Subgroup numbers are variant-specific

The same integer does not mean the same fragment in different variants. In the
original UNIFAC set, subgroup 118 is the sulfone fragment `(CH2)2SU`; in the
PSRK set, subgroup 118 is methane `CH4`. Passing a number from one variant's
table into another variant's constructor therefore builds a model for a
different chemistry, and the program will not complain: the lookup falls back
to whatever the parameter table happens to hold, and the results are silently
wrong.

Two practical consequences:

- Always take subgroup numbers from the table of the variant you are
  constructing. Each variant page in this section has its own table.
- If a result looks absurd (a coefficient of order 10^-700, or a coefficient
  that jumps by orders of magnitude between two similar mixtures), suspect a
  misidentified or out-of-table subgroup before suspecting the physics.

Molecular size and area parameters are sums of subgroup contributions,
\(r_i=\sum_k\nu_k^i R_k\) and \(q_i=\sum_k\nu_k^i Q_k\), where
\(\nu_k^i\) is the count of subgroup `k` in component `i`. `R_k` and `Q_k` are
dimensionless relative volumes and areas, not physical SI measurements.

The complete example below can be used without reading the derivative reference.
The reference explains how the implementation obtains thermal and compositional
derivatives consistently from the same excess-energy function.

## Combinatorial term derivatives

Combinatorial term has two parameters \(z\) and \(d\) that could be modified.
The \(z\) parameter is always set to 10, and the \(d\) parameter is set to 1
for the Classic Liquid-Vapor UNIFAC model and the PSRK-UNIFAC model, and 3/4
for the Dortmund-UNIFAC model. Both could be changed (they are attributes of
the UNIFAC objects (`z` and `d` respectively)), but changing them creates a
custom model: it does not transfer the published validation of another variant.

Calculate the UNIFAC combinatorial term of reduced Gibbs excess energy.
The subroutine uses the Flory-Huggins and Staverman-Guggenheim
combinatory terms as follows:

### Flory-Huggins

$$
   G^{E,FH} =
   RT \left(\sum_i^{NC} n_i \, \text{ln} \, r_i^d
   - n \, \text{ln} \, \sum_j^{NC} n_j r_j^d
   + n \, \text{ln} \, n \right)
$$

$$
   \frac{dG^{E,FH}}{dn_i} =
   RT \left(\text{ln} \, r_i^d - \text{ln} \, \sum_j^{NC} n_j r_j^d
   + \text{ln} \, n + 1 - \frac{n r_i^d}{\displaystyle
   \sum_j^{NC} n_j r_j^d} \right)
$$

$$
   \frac{d^2G^{E,FH}}{dn_i dn_j} =
   RT \left(- \frac{r_i^d + r_j^d}{\displaystyle \sum_l^{NC} n_l r_l^d}
   + \frac{1}{n} + \frac{n r_i^d r_j^d}{\displaystyle \left(\sum_l^{NC}
   n_l r_l^d \right)^2} \right)
$$

### Staverman-Guggenheim

$$
   \frac{G^{E,SG}}{RT} =
   \frac{z}{2} \sum_i^{NC} n_i q_i
   \left(\text{ln} \frac{q_i}{r_i}
   - \text{ln} \, \sum_j^{NC} n_j q_j
   + \text{ln} \, \sum_j^{NC} n_j r_j \right)
$$

$$
   \frac{1}{RT}\frac{dG^{E,SG}}{dn_i} =
   \frac{z}{2} q_i \left(
   - \text{ln} \, \left(
   \frac{r_i \sum_j^{NC} n_j q_j}{\displaystyle q_i \sum_j^{NC}
   n_j r_j} \right) - 1 + \frac{\displaystyle r_i \sum_j^{NC} n_j
   q_j}{\displaystyle q_i \sum_j^{NC} n_j r_j} \right)
$$

$$
   \frac{1}{RT}\frac{d^2G^{E,SG}}{dn_i dn_j} =
   \frac{z}{2} \left(- \frac{q_i q_j}{\displaystyle \sum_l^{NC} n_lq_l}
   + \frac{q_i r_j + q_j r_i}{\displaystyle \sum_l^{NC} n_l r_l}
   - \frac{\displaystyle r_i r_j \sum_l^{NC} n_l q_l}
   {\left(\displaystyle \sum_l^{NC} n_l r_l \right)^2} \right)
$$

### Fredenslund et al. (UNIFAC)

The complete combinatorial term is given by the sum of the Flory-Huggins and
Staverman-Guggenheim terms:

$$
   \frac{G^{E,\text{UNIFAC}}}{RT} =
   \frac{G^{E,FH}}{RT} + \frac{G^{E,SG}}{RT}
$$

## Residual terms derivatives
Evaluate the UNIFAC residual term. Check the temperature function \(E_{jk}\).
This temperature function is different for each UNIFAC model. The residual
Gibbs excess energy and its derivatives are evaluated as:

$$
 \frac{G^{E,R}}{RT} = - \sum_i^{NC} n_i \sum_k^{NG} v_k^i Q_k
 (\Lambda_k - \Lambda_k^i)
$$

With:

$$
 \Lambda_k = \text{ln} \, \sum_{j}^{NG} \Theta_j E_{jk}
$$

$$
 \Lambda_k^i = \text{ln} \, \sum_{j}^{NG} \Theta_j^i E_{jk}
$$

$$
 E_{jk} = \text{exp} \left(- \frac{U_{jk}}{RT} \right)
$$

$$
 \Theta_j = \frac{Q_j \displaystyle \sum_{l}^{NC} n_l v_j^l}
 {\displaystyle \sum_{l}^{NC} n_l \sum_{m}^{NG} v_m^l Q_m}
$$

$$
 \Theta_j^i = \frac{Q_j v_j^i}{\displaystyle \sum_k^{NG} v_k^i Q_k}
$$

In the UNIFAC model, the \(\Theta_j^i \) values are calculated assuming
that the molecule "i" is pure, hence only the subgroups of the molecule
"i" must be considered for the calculation. On the other hand, for the
\(\Theta_j \) values, all the system's subgroups are considered.

### The compositional derivatives:

$$
 \frac{1}{R T} \frac{\partial G^{E,R}}{\partial n_\alpha} =
 - \sum_k^{\mathrm{NG}} v_k^\alpha Q_k \left(\Lambda_k -
 \Lambda_k^\alpha \right) - \sum_i^{\mathrm{NC}} n_i
 \sum_k^{\mathrm{NG}} v_k^i Q_k
 \frac{\partial \Lambda_k}{\partial n_\alpha}
$$

$$
 \frac{1}{R T} \frac{\partial^2 G^{E,R}}{\partial n_
 \alpha \partial n_\beta} = -\sum_k^{\mathrm{NG}} Q_k \left(v_k^\alpha
 \frac{\partial \Lambda_k}{\partial n_\beta} + v_k^\beta
 \frac{\partial \Lambda_k}{\partial n_\alpha}\right)
 - \sum_k^{\mathrm{NG}} \left(\sum_i^{\mathrm{NC}} n_i v_k^i\right) Q_k
 \frac{\partial^2 \Lambda_k}{\partial n_\alpha \partial n_\beta}
$$

With:

$$
 \frac{\partial \Lambda_k}{\partial n_\alpha}
 = \frac{\sum_j^{\mathrm{NG}} v_j^\alpha Q_j E_{j k}}
 {\sum_l^{\mathrm{NC}} n_l \sum_j^{\mathrm{NG}} v_j^l Q_j
 E_{j k}} - \frac{\sum_m^{\mathrm{NG}} v_m^\alpha Q_m}
 {\sum_l^{\mathrm{NC}} n_l \sum_m^{\mathrm{NG}} v_m^l Q_m}
$$

$$
 \frac{\partial^2 \Lambda_k}{\partial n_\alpha \partial n_\beta}
 = - \frac{\left(\sum_j^{\mathrm{NG}} v_j^\alpha Q_j E_{j k}\right)
 \left(\sum_j^{\mathrm{NG}} v_j^\beta Q_j E_{j k}\right)}
 {\left(\sum_l^{\mathrm{NC}} n_l \sum_j^{\mathrm{NG}} v_j^l Q_j
 E_{j k}\right)^2} + \frac{\left(\sum_m^{\mathrm{NG}} v_m^\alpha
 Q_m\right)\left(\sum_m^{\mathrm{NG}} v_m^\beta Q_m\right)}
 {\left(\sum_l^{\mathrm{NC}} n_l
 \sum_m^{\mathrm{NG}} v_m^l Q_m\right)^2}
$$

### The temperature derivatives:

$$
 \frac{\partial\left(\frac{G^{E, R}}{R T}\right)}{\partial T} =
 -\sum_i^{\mathrm{NC}} n_i \sum_k^{\mathrm{NG}} v_k^i Q_k
 \left(\frac{\partial \Lambda_k}{\partial T}
 -\frac{\partial \Lambda_k^i}{\partial T}\right)
$$

$$
 \frac{\partial^2\left(\frac{G^{E,R}}{R T}\right)}{\partial T^2} =
 -\sum_i^{\mathrm{NC}} n_i \sum_k^{\mathrm{NG}} v_k^i Q_k
 \left(\frac{\partial^2 \Lambda_k}{\partial T^2} -
 \frac{\partial^2 \Lambda_k^i}{\partial T^2}\right)
$$

With:

$$
 \frac{\partial \Lambda_k}{\partial T} =
 \frac{\sum_{j}^{NG} \Theta_j \frac{d E_{jk}}{dT}}
 {\sum_{j}^{NG} \Theta_j E_{jk}}
$$

$$
 \frac{\partial \Lambda_k^i}{\partial T} =
 \frac{\sum_{j}^{NG} \Theta_j^i \frac{d E_{jk}}{dT}}
 {\sum_{j}^{NG} \Theta_j^i E_{jk}}
$$

$$
 \frac{\partial^2 \Lambda_k}{\partial T^2} =
 \frac{\sum_{j}^{NG} \Theta_j \frac{d^2 E_{jk}}{dT^2}}
 {\sum_{j}^{NG} \Theta_j E_{jk}}
 - \left(\frac{\partial \Lambda_k}{\partial T} \right)^2
$$

$$
 \frac{\partial^2 \Lambda_k^i}{\partial T^2} =
 \frac{\sum_{j}^{NG} \Theta_j^i \frac{d^2 E_{jk}}{dT^2}}
 {\sum_{j}^{NG} \Theta_j^i E_{jk}}
 - \left(\frac{\partial \Lambda_k^i}{\partial T} \right)^2
$$

### Temperature-compositional cross derivative:

$$
 \frac{\partial \left(\frac{G^{E, R}}{R T} \right)}
 {\partial n_\alpha \partial T}=
 -\sum_k^{\mathrm{NG}} v_k^\alpha Q_k \left(\frac{\partial \Lambda_k}
 {\partial T} - \frac{\partial \Lambda_k^\alpha}{\partial T}\right)
 -\sum_k^{\mathrm{NG}} \left(\sum_i^{\mathrm{NC}} n_i v_k^i \right)
 Q_k \frac{\partial^2 \Lambda_k}{\partial n_\alpha \partial T}
$$

With:

$$
  \frac{\partial^2 \Lambda_k}{\partial n_\alpha \partial T} =
 \frac{\sum_j^{\mathrm{NG}} v_j^\alpha Q_j \frac{\partial
 \tilde{E}_{j k}}{\partial T}}{\sum_l^{\mathrm{NC}} n_l
 \sum_j^{\mathrm{NG}} v_j^l Q_j \tilde{E}_{j k}} -
 \frac{\left(\sum_j^{\mathrm{NG}} v_j^\alpha Q_j \tilde{E}_{j k}\right)
 \left(\sum_l^{\mathrm{NC}} n_l \sum_j^{\mathrm{NG}} v_j^l Q_j
 \frac{\partial \tilde{E}_{j k}}{\partial T}\right)}
 {\left(\sum_l^{\mathrm{NC}} n_l
 \sum_j^{\mathrm{NG}} v_j^l Q_j \tilde{E}_{j k}\right)^2}
$$

## Examples

On each specific UNIFAC model section we will give more examples of each
specific model. Here we leave a simple example for the classic Liquid-Vapor
UNIFAC model.

### Calculating activity coefficients
We can instantiate a [[UNIFAC]] model with a mixture ethanol-water and evaluate
the logarithm of activity coefficients of the model for a 0.5 mole fraction of
each, and a temperature of 298.0 K.

```fortran
program unifac_activity
use yaeos, only: pr, Groups, UNIFAC, setup_unifac
implicit none

! Variables declarations
type(UNIFAC) :: model
type(Groups) :: molecules(2)
real(pr) :: ln_gammas(2)
real(pr) :: n(2)
real(pr) :: T

! Variables instances
! Ethanol definition [CH3, CH2, OH]
molecules(1)%groups_ids = [1, 2, 14] ! Subgroups ids
molecules(1)%number_of_groups = [1, 1, 1] ! Subgroups occurrences

! Water definition [H2O]
molecules(2)%groups_ids = [16]
molecules(2)%number_of_groups = [1]

! Model setup
model = setup_unifac(molecules)

! Calculate ln_gammas
n = [0.5_pr, 0.5_pr]
T = 298.0_pr
call model%ln_activity_coefficient(n, T, ln_gammas)

print *, "ln(gamma): ethanol, water = ", ln_gammas
print *, "gamma = ", exp(ln_gammas)
print *, "activity = ", n/sum(n)*exp(ln_gammas)
end program unifac_activity
```

This is a complete program. Each output array is ordered ethanol, then water.
`ln_gammas` are natural logarithms, `exp(ln_gammas)` are coefficients, and the
final line multiplies by the mole fractions to obtain activities. Positive
logarithms mean coefficients greater than 1; neither coefficients nor activities
should be interpreted as percentages.

### Worked example: the same mixture through two variants

Because the three variants share the same Fortran type, you can compare them
with one extra program. This is a good way to see what the "variant" choice
actually changes.

```fortran
program unifac_variants
use yaeos, only: pr, Groups, UNIFAC
use yaeos, only: setup_unifac, setup_dortmund, setup_psrk
implicit none

type(UNIFAC) :: model_lv
type(UNIFAC) :: model_do
type(UNIFAC) :: model_ps
type(Groups) :: molecules(2)
real(pr) :: ln_lv(2)
real(pr) :: ln_do(2)
real(pr) :: ln_ps(2)
real(pr) :: n(2)
real(pr) :: T

molecules(1)%groups_ids = [1, 2, 14]
molecules(1)%number_of_groups = [1, 1, 1]
molecules(2)%groups_ids = [16]
molecules(2)%number_of_groups = [1]

n = [0.5_pr, 0.5_pr]
T = 298.0_pr

model_lv = setup_unifac(molecules)
model_do = setup_dortmund(molecules)
model_ps = setup_psrk(molecules)

call model_lv%ln_activity_coefficient(n, T, ln_lv)
call model_do%ln_activity_coefficient(n, T, ln_do)
call model_ps%ln_activity_coefficient(n, T, ln_ps)

print *, "original UNIFAC ln(gamma) = ", ln_lv
print *, "Dortmund     ln(gamma) = ", ln_do
print *, "PSRK         ln(gamma) = ", ln_ps
end program unifac_variants
```

The Dortmund answer differs from the other two. That is expected: each variant
has its own fitted group interactions and its own geometry, so in general the
numbers are independent estimates rather than one number computed several ways.
For many ordinary mixtures the PSRK parameter set reuses the same size
parameters and `a` values as the original model for the common main groups;
for this ethanol-water example the original and PSRK results coincide, while
the Dortmund set, with its modified size term and fitted interactions, differs.
A visible difference is a reminder that these are predictions from different
data; it is not a bug.

### Worked example: a liquid-liquid split

Activity models are usually used through the phase-equilibrium routines rather
than directly. The program below uses UNIFAC to find a two-liquid split of a
ternary feed, which is the classic use case for a group-contribution model:
screening which components separate without binary-specific fits.

```fortran
program unifac_lle
use yaeos, only: pr, Groups, UNIFAC, setup_unifac
use yaeos, only: EquilibriumState, flash, min_tpd
implicit none

integer, parameter :: nc = 3

type(UNIFAC) :: model
type(Groups) :: molecules(nc)
type(EquilibriumState) :: result
real(pr) :: n(nc)
real(pr) :: w(nc)
real(pr) :: T
real(pr) :: P
real(pr) :: mintpd
integer :: iters

molecules(1)%groups_ids = [1, 42]
molecules(1)%number_of_groups = [1, 1]

molecules(2)%groups_ids = [1, 2]
molecules(2)%number_of_groups = [2, 4]

molecules(3)%groups_ids = [16]
molecules(3)%number_of_groups = [1]

model = setup_unifac(molecules)

n = [0.2_pr, 0.7_pr, 0.1_pr]
T = 250.0_pr
P = 1.0_pr

call min_tpd(model, n, P=P, T=T, mintpd=mintpd, w=w)
print *, "mintpd = ", mintpd

if (mintpd < 0.0_pr) then
   result = flash(model, n, T, k0=w/n, iters=iters)
   print *, "kind = ", result%kind
   print *, "x = ", result%x
   print *, "y = ", result%y
   print *, "beta = ", result%beta
else
   print *, "Feed is stable in this search."
end if
end program unifac_lle
```

`min_tpd` looks for an instability of the feed; a negative value means the
mixture wants to split. `flash` then solves the isofugacity and mass-balance
equations at the given temperature, starting from the trial composition. For an
activity model it requires the initial K-factors `k0` and cannot accept a
pressure or volume specification. Check `result%kind` before using `x`, `y`,
and `beta`; the pressure and phase volumes are not set for a `GeModel` result.

The [flash page](../../phase_equilibrium/flash.html) describes this workflow in
more detail, including how to initialize from a previous nearby calculation
when building tables of phase equilibria.

## Checks, limitations, and next steps

A few quick checks catch most mistakes before they become wrong results.

**Check the amount scaling.** Multiply every amount by a common factor and
confirm that the activity coefficients do not change, while the extensive
`Ge`, `He`, and `Se` scale with that factor. A model that fails this test has
a normalization problem.

**Check the pure-component limit.** Approach a pure composition and confirm
that the coefficient of the *present* component tends to 1. This does not
require a trace solute's coefficient to approach 1: a solute at infinite
dilution has its own limiting coefficient. The pure-group reference subtraction
in the residual term enforces this standard-state convention, and it matters
even for molecules made of several groups.

**Check the decomposition before the physics.** When results look wrong, look
at the group assignment first. A wrong OH subgroup, an omitted ring-specific
fragment, or a double-counted compound fragment changes both the geometric and
the energetic contribution at once, so there is no way to detect it from the
output alone.

**Know the limits.** Prediction quality depends on the fitted coverage of the
groups involved and on the temperature range of that fit. More temperature
coefficients do not guarantee reliable extrapolation. The standard variants are
not general electrolyte or speciation models, and neither liquid density nor
vapor properties follow from an isolated UNIFAC object.

Continue with the variant's subgroup tables and examples, then use the
[shared excess-property methods](../index.html) for thermal properties and
derivatives. Compare with experimental data before drawing process-design or
phase-split conclusions. A trial liquid's coefficients alone do not establish
its stability; see [phase equilibrium](../../phase_equilibrium/index.html).

# References
1. [Dortmund Data Bank Software & Separation Technology](https://www.ddbst.com/published-parameters-unifac.html)
2. Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975). Group‐contribution
estimation of activity coefficients in nonideal liquid mixtures. AIChE Journal,
21(6), 1086–1099.
[https://doi.org/10.1002/aic.690210607](https://doi.org/10.1002/aic.690210607)
3. Skjold-Jorgensen, S., Kolbe, B., Gmehling, J., & Rasmussen, P. (1979).
Vapor-Liquid Equilibria by UNIFAC Group Contribution. Revision and Extension.
Industrial & Engineering Chemistry Process Design and Development, 18(4),
714–722.
[https://doi.org/10.1021/i260072a024](https://doi.org/10.1021/i260072a024)
4. Gmehling, J., Rasmussen, P., & Fredenslund, A. (1982). Vapor-liquid
equilibriums by UNIFAC group contribution. Revision and extension. 2.
Industrial & Engineering Chemistry Process Design and Development, 21(1),
118–127.
[https://doi.org/10.1021/i200016a021](https://doi.org/10.1021/i200016a021)
5. Macedo, E. A., Weidlich, U., Gmehling, J., & Rasmussen, P. (1983).
Vapor-liquid equilibriums by UNIFAC group contribution. Revision and extension. 3.
Industrial & Engineering Chemistry Process Design and Development, 22(4),
676–678.
[https://doi.org/10.1021/i200023a023](https://doi.org/10.1021/i200023a023)
6. Tiegs, D., Rasmussen, P., Gmehling, J., & Fredenslund, A. (1987).
Vapor-liquid equilibria by UNIFAC group contribution. 4. Revision and
extension. Industrial & Engineering Chemistry Research, 26(1), 159–161.
[https://doi.org/10.1021/ie00061a030](https://doi.org/10.1021/ie00061a030)
7. Hansen, H. K., Rasmussen, P., Fredenslund, A., Schiller, M., & Gmehling, J.
(1991). Vapor-liquid equilibria by UNIFAC group contribution. 5. Revision and
extension. Industrial & Engineering Chemistry Research, 30 (10), 2352–2355.
[https://doi.org/10.1021/ie00058a017](https://doi.org/10.1021/ie00058a017)
8. Wittig, R., Lohmann, J., & Gmehling, J. (2003). Vapor−Liquid Equilibria by
UNIFAC Group Contribution. 6. Revision and Extension. Industrial & Engineering
Chemistry Research, 42(1), 183–188.
[https://doi.org/10.1021/ie020506l](https://doi.org/10.1021/ie020506l)
9. [SINTEF - Thermopack](https://github.com/termotools/thermopack)
