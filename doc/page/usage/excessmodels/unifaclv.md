---
title: UNIFAC-LV
---

[TOC]

# UNIFAC

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

Each substance of a mixture modeled with [[UNIFAC]] must be represented by a
list a functional groups and other list with the ocurrence of each functional
group on the substance. The list of the functional groups culd be accesed on
the DDBST web page:
[https://www.ddbst.com/published-parameters-unifac.html](https://www.ddbst.com/published-parameters-unifac.html)

## Examples
### Calculating activity coefficients
We can instantiate a [[UNIFAC]] model with a mixture ethanol-water and evaluate
the logarithm of activity coefficients of the model for a 0.5 mole fraction of
each, and a temperature of 298.0 K.

```fortran
use yaeos__constants, only: pr
use yaeos__models_ge_group_contribution_unifac, only: Groups, UNIFAC, setup_unifac

! Variables declarations
type(UNIFAC) :: model
type(Groups) :: molecules(2)
real(pr) :: ln_gammas(2)

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
call model%ln_activity_coefficient([0.5_pr, 0.5_pr], 298.0_pr, ln_gammas)

print *, ln_gammas
```

You will obtain:

```
>>> 0.18534142000449058    0.40331395945417559
```

# References
1. [Dortmund Data Bank Software & Separation Technology](https://www.ddbst
.com/published-parameters-unifac.html)
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
Vapor-liquid equilibriums by UNIFAC group contribution. Revision and extension.
3. Industrial & Engineering Chemistry Process Design and Development, 22(4),
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
9. [SINTEF - Thermopack](https://github.com/thermotools/thermopack)

