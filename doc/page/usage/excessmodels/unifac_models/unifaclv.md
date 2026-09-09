---
title: UNIFAC-LV
---

[TOC]

# Liquid-Vapor UNIFAC

## Model description

This is the original and classic Liquid-Vapor UNIFAC model. In this model, the
parameters are defined as:
$$
    z = 10
$$

$$
    d = 1
$$

The temperature function \(E_{jk}\) is defined with a single temperature 
constant \(a_{jk}\) coefficient as follows:

$$
    E_{jk} = \text{exp} \left(- \frac{U_{jk}}{RT} \right) = \text{exp} \left(- \frac{a_{jk}}{T} \right)
$$

## Subgroups list

The list of the functional groups and its interaction parameters could be 
accessed on the DDBST web page:
[https://www.ddbst.com/published-parameters-unifac.html](https://www.ddbst.com/published-parameters-unifac.html)

We reproduce here the list of functional groups. To instantiate a UNIFAC model
you must define which functional groups are used in a molecule by the Subgroup
Number column.

|   Subgroup number | Subgroup Name   | Maingroup    |      R |      Q |
|:-----------------:|:---------------:|:------------:|:------:|:------:|
|                 1 | CH3             | [1]CH2       | 0.9011 | 0.848  |
|                 2 | CH2             | [1]CH2       | 0.6744 | 0.54   |
|                 3 | CH              | [1]CH2       | 0.4469 | 0.228  |
|                 4 | C               | [1]CH2       | 0.2195 | 0      |
|                 5 | CH2=CH          | [2]C=C       | 1.3454 | 1.176  |
|                 6 | CH=CH           | [2]C=C       | 1.1167 | 0.867  |
|                 7 | CH2=C           | [2]C=C       | 1.1173 | 0.988  |
|                 8 | CH=C            | [2]C=C       | 0.8886 | 0.676  |
|                 9 | ACH             | [3]ACH       | 0.5313 | 0.4    |
|                10 | AC              | [3]ACH       | 0.3652 | 0.12   |
|                11 | ACCH3           | [4]ACCH2     | 1.2663 | 0.968  |
|                12 | ACCH2           | [4]ACCH2     | 1.0396 | 0.66   |
|                13 | ACCH            | [4]ACCH2     | 0.8121 | 0.348  |
|                14 | OH              | [5]OH        | 1      | 1.2    |
|                15 | CH3OH           | [6]CH3OH     | 1.4311 | 1.432  |
|                16 | H2O             | [7]H2O       | 0.92   | 1.4    |
|                17 | ACOH            | [8]ACOH      | 0.8952 | 0.68   |
|                18 | CH3CO           | [9]CH2CO     | 1.6724 | 1.488  |
|                19 | CH2CO           | [9]CH2CO     | 1.4457 | 1.18   |
|                20 | HCO             | [10]HCO      | 0.998  | 0.948  |
|                21 | CH3COO          | [11]CCOO     | 1.9031 | 1.728  |
|                22 | CH2COO          | [11]CCOO     | 1.6764 | 1.42   |
|                23 | HCOO            | [12]HCOO     | 1.242  | 1.188  |
|                24 | CH3O            | [13]CH2O     | 1.145  | 1.088  |
|                25 | CH2O            | [13]CH2O     | 0.9183 | 0.78   |
|                26 | CHO             | [13]CH2O     | 0.6908 | 0.468  |
|                27 | THF             | [13]CH2O     | 0.9183 | 1.1    |
|                28 | CH3NH2          | [14]CNH2     | 1.5959 | 1.544  |
|                29 | CH2NH2          | [14]CNH2     | 1.3692 | 1.236  |
|                30 | CHNH2           | [14]CNH2     | 1.1417 | 0.924  |
|                31 | CH3NH           | [15]CNH      | 1.4337 | 1.244  |
|                32 | CH2NH           | [15]CNH      | 1.207  | 0.936  |
|                33 | CHNH            | [15]CNH      | 0.9795 | 0.624  |
|                34 | CH3N            | [16] (C)3N   | 1.1865 | 0.94   |
|                35 | CH2N            | [16] (C)3N   | 0.9597 | 0.632  |
|                36 | ACNH2           | [17]ACNH2    | 1.06   | 0.816  |
|                37 | C5H5N           | [18]PYRIDINE | 2.9993 | 2.113  |
|                38 | C5H4N           | [18]PYRIDINE | 2.8332 | 1.833  |
|                39 | C5H3N           | [18]PYRIDINE | 2.667  | 1.553  |
|                40 | CH3CN           | [19]CCN      | 1.8701 | 1.724  |
|                41 | CH2CN           | [19]CCN      | 1.6434 | 1.416  |
|                42 | COOH            | [20]COOH     | 1.3013 | 1.224  |
|                43 | HCOOH           | [20]COOH     | 1.528  | 1.532  |
|                44 | CH2CL           | [21]CCL      | 1.4654 | 1.264  |
|                45 | CHCL            | [21]CCL      | 1.238  | 0.952  |
|                46 | CCL             | [21]CCL      | 1.0106 | 0.724  |
|                47 | CH2CL2          | [22]CCL2     | 2.2564 | 1.988  |
|                48 | CHCL2           | [22]CCL2     | 2.0606 | 1.684  |
|                49 | CCL2            | [22]CCL2     | 1.8016 | 1.448  |
|                50 | CHCL3           | [23]CCL3     | 2.87   | 2.41   |
|                51 | CCL3            | [23]CCL3     | 2.6401 | 2.184  |
|                52 | CCL4            | [24]CCL4     | 3.39   | 2.91   |
|                53 | ACCL            | [25]ACCL     | 1.1562 | 0.844  |
|                54 | CH3NO2          | [26]CNO2     | 2.0086 | 1.868  |
|                55 | CH2NO2          | [26]CNO2     | 1.7818 | 1.56   |
|                56 | CHNO2           | [26]CNO2     | 1.5544 | 1.248  |
|                57 | ACNO2           | [27]ACNO2    | 1.4199 | 1.104  |
|                58 | CS2             | [28]CS2      | 2.057  | 1.65   |
|                59 | CH3SH           | [29]CH3SH    | 1.877  | 1.676  |
|                60 | CH2SH           | [29]CH3SH    | 1.651  | 1.368  |
|                61 | FURFURAL        | [30]FURFURAL | 3.168  | 2.484  |
|                62 | DOH             | [31]DOH      | 2.4088 | 2.248  |
|                63 | I               | [32]I        | 1.264  | 0.992  |
|                64 | BR              | [33]BR       | 0.9492 | 0.832  |
|                65 | CH=-C           | [34]C=-C     | 1.292  | 1.088  |
|                66 | C=-C            | [34]C=-C     | 1.0613 | 0.784  |
|                67 | DMSO            | [35]DMSO     | 2.8266 | 2.472  |
|                68 | ACRY            | [36]ACRY     | 2.3144 | 2.052  |
|                69 | CL-(C=C)        | [37]CLCC     | 0.791  | 0.724  |
|                70 | C=C             | [2]C=C       | 0.6605 | 0.485  |
|                71 | ACF             | [38]ACF      | 0.6948 | 0.524  |
|                72 | DMF             | [39]DMF      | 3.0856 | 2.736  |
|                73 | HCON(CH2)2      | [39]DMF      | 2.6322 | 2.12   |
|                74 | CF3             | [40]CF2      | 1.406  | 1.38   |
|                75 | CF2             | [40]CF2      | 1.0105 | 0.92   |
|                76 | CF              | [40]CF2      | 0.615  | 0.46   |
|                77 | COO             | [41]COO      | 1.38   | 1.2    |
|                78 | SIH3            | [42]SIH2     | 1.6035 | 1.2632 |
|                79 | SIH2            | [42]SIH2     | 1.4443 | 1.0063 |
|                80 | SIH             | [42]SIH2     | 1.2853 | 0.7494 |
|                81 | SI              | [42]SIH2     | 1.047  | 0.4099 |
|                82 | SIH2O           | [43]SIO      | 1.4838 | 1.0621 |
|                83 | SIHO            | [43]SIO      | 1.303  | 0.7639 |
|                84 | SIO             | [43]SIO      | 1.1044 | 0.4657 |
|                85 | NMP             | [44]NMP      | 3.981  | 3.2    |
|                86 | CCL3F           | [45]CCLF     | 3.0356 | 2.644  |
|                87 | CCL2F           | [45]CCLF     | 2.2287 | 1.916  |
|                88 | HCCL2F          | [45]CCLF     | 2.406  | 2.116  |
|                89 | HCCLF           | [45]CCLF     | 1.6493 | 1.416  |
|                90 | CCLF2           | [45]CCLF     | 1.8174 | 1.648  |
|                91 | HCCLF2          | [45]CCLF     | 1.967  | 1.828  |
|                92 | CCLF3           | [45]CCLF     | 2.1721 | 2.1    |
|                93 | CCL2F2          | [45]CCLF     | 2.6243 | 2.376  |
|                94 | AMH2            | [46]CON(AM)  | 1.4515 | 1.248  |
|                95 | AMHCH3          | [46]CON(AM)  | 2.1905 | 1.796  |
|                96 | AMHCH2          | [46]CON(AM)  | 1.9637 | 1.488  |
|                97 | AM(CH3)2        | [46]CON(AM)  | 2.8589 | 2.428  |
|                98 | AMCH3CH2        | [46]CON(AM)  | 2.6322 | 2.12   |
|                99 | AM(CH2)2        | [46]CON(AM)  | 2.4054 | 1.812  |
|               100 | C2H5O2          | [47]OCCOH    | 2.1226 | 1.904  |
|               101 | C2H4O2          | [47]OCCOH    | 1.8952 | 1.592  |
|               102 | CH3S            | [48]CH2S     | 1.613  | 1.368  |
|               103 | CH2S            | [48]CH2S     | 1.3863 | 1.06   |
|               104 | CHS             | [48]CH2S     | 1.1589 | 0.748  |
|               105 | MORPH           | [49]MORPH    | 3.474  | 2.796  |
|               106 | C4H4S           | [50]THIOPHEN | 2.8569 | 2.14   |
|               107 | C4H3S           | [50]THIOPHEN | 2.6908 | 1.86   |
|               108 | C4H2S           | [50]THIOPHEN | 2.5247 | 1.58   |
|               109 | NCO             | [51]NCO      | 1.0567 | 0.732  |
|               118 | (CH2)2SU        | [55]SULFONES | 2.6869 | 2.12   |
|               119 | CH2CHSU         | [55]SULFONES | 2.4595 | 1.808  |
|               178 | IMIDAZOL        | [84]IMIDAZOL | 2.026  | 0.868  |
|               179 | BTI             | [85]BTI      | 5.774  | 4.932  |


## Using `ugropy` to retrieve UNIFAC subgroups

There is the possibility of using another library of our group
[`ugropy`](https://github.com/ipqa-research/ugropy) to retrieve the UNIFAC 
subgroups and not suffer the pain of typing the subgroup numbers and
parameters by hand. The next Python snippet shows how you can use it.

```python
from ugropy import unifac, writers


names = ["water", "toluene", "acetone"]
groups = [unifac.get_groups(n).subgroups for n in names]

fortran_code = writers.to_yaeos(groups, unifac)

print(fortran_code)
```

And you will obtain:

```fortran
use yaeos__models_ge_group_contribution_unifac, only: Groups

type(Groups) :: molecules(3)

molecules(1)%groups_ids = [16]
molecules(1)%number_of_groups = [1]

molecules(2)%groups_ids = [9, 11]
molecules(2)%number_of_groups = [5, 1]

molecules(3)%groups_ids = [1, 18]
molecules(3)%number_of_groups = [1, 1]
```


## Examples

Here is an example of a fully instantiated UNIFACLV model. Please check the
`Gibbs Excess Models` section in the user documentation to learn all the things
you can do with this model.

Notice that here we are using the [[setup_unifac]] function to instantiate the
model.

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

## References

1. Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975). Group‐contribution
estimation of activity coefficients in nonideal liquid mixtures. AIChE Journal,
21(6), 1086–1099.
[https://doi.org/10.1002/aic.690210607](https://doi.org/10.1002/aic.690210607)
2. Skjold-Jorgensen, S., Kolbe, B., Gmehling, J., & Rasmussen, P. (1979).
Vapor-Liquid Equilibria by UNIFAC Group Contribution. Revision and Extension.
Industrial & Engineering Chemistry Process Design and Development, 18(4),
714–722.
[https://doi.org/10.1021/i260072a024](https://doi.org/10.1021/i260072a024)
3. Gmehling, J., Rasmussen, P., & Fredenslund, A. (1982). Vapor-liquid
equilibriums by UNIFAC group contribution. Revision and extension. 2.
Industrial & Engineering Chemistry Process Design and Development, 21(1),
118–127.
[https://doi.org/10.1021/i200016a021](https://doi.org/10.1021/i200016a021)
4. Macedo, E. A., Weidlich, U., Gmehling, J., & Rasmussen, P. (1983).
Vapor-liquid equilibriums by UNIFAC group contribution. Revision and extension.
5. Industrial & Engineering Chemistry Process Design and Development, 22(4),
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