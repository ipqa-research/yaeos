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

## Combinatorial term derivatives

Combinatorial term has two parameters \(z\) and \(d\) that could be modified.
The \(z\) parameter is always set to 10, and the \(d\) parameter is set to 1
for the Classic Liquid-Vapor UNIFAC model and the PSRK-UNIFAC model, and 3/4
for the Dortmund-UNIFAC model. Both could be changed (they are attributes of
the UNIFAC objects).

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
$$
   \frac{G^{E,\text{UNIFAC}}}{RT} =
   \frac{G^{E,FH}}{RT} + \frac{G^{E,SG}}{RT}
$$

## Residual terms derivatives
Evaluate the UNIFAC residual term. The residual Gibbs excess energy
and its derivatives are evaluated as:

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
 {\displaystyle \sum_{k}^{NC} n_k \sum_{m}^{NG} v_m^l Q_m}
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
9. [SINTEF - Thermopack](https://github.com/termotools/thermopack)

