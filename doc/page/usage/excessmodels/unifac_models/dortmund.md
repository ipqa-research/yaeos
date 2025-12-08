---
title: UNIFAC-Dortmund
---

[TOC]

# Dortmund modified UNIFAC model

## Model description

This is the Dortmund modified UNIFAC model. In this model, the parameters are
defined as:

$$
    z = \frac{3}{4}
$$

$$
    d = 1
$$

The temperature function \(E_{jk}\) is defined with a quadratic temperature
function as follows:

$$
    E_{jk} = \text{exp} \left(- \frac{U_{jk}}{RT} \right) = 
    \text{exp} \left(- \frac{a_{jk} + b_{jk} T + c_{jk} T^2}{T} \right)
$$

## Subgroups list

The list of the functional groups and its interaction parameters could be 
accessed on the DDBST web page:
[https://www.ddbst.com/PublishedParametersUNIFACDO.html](https://www.ddbst.com/PublishedParametersUNIFACDO.html)

We reproduce here the list of functional groups. To instantiate a
UNIFAC-Dortmund model you must define which functional groups are used in a
molecule by the Subgroup Number column.


|   Subgroup number | Subgroup Name   |   Main Group No. | Main Group Name   |      R |      Q |
|:-----------------:|:---------------:|:----------------:|:-----------------:|:------:|:------:|
|                 1 | CH3             |                1 | CH2               | 0.6325 | 1.0608 |
|                 2 | CH2             |                1 | CH2               | 0.6325 | 0.7081 |
|                 3 | CH              |                1 | CH2               | 0.6325 | 0.3554 |
|                 4 | C               |                1 | CH2               | 0.6325 | 0      |
|                 5 | CH2=CH          |                2 | C=C               | 1.2832 | 1.6016 |
|                 6 | CH=CH           |                2 | C=C               | 1.2832 | 1.2489 |
|                 7 | CH2=C           |                2 | C=C               | 1.2832 | 1.2489 |
|                 8 | CH=C            |                2 | C=C               | 1.2832 | 0.8962 |
|                 9 | ACH             |                3 | ACH               | 0.3763 | 0.4321 |
|                10 | AC              |                3 | ACH               | 0.3763 | 0.2113 |
|                11 | ACCH3           |                4 | ACCH2             | 0.91   | 0.949  |
|                12 | ACCH2           |                4 | ACCH2             | 0.91   | 0.7962 |
|                13 | ACCH            |                4 | ACCH2             | 0.91   | 0.3769 |
|                14 | OH (P)          |                5 | OH                | 1.2302 | 0.8927 |
|                15 | CH3OH           |                6 | CH3OH             | 0.8585 | 0.9938 |
|                16 | H2O             |                7 | H2O               | 1.7334 | 2.4561 |
|                17 | ACOH            |                8 | ACOH              | 1.08   | 0.975  |
|                18 | CH3CO           |                9 | CH2CO             | 1.7048 | 1.67   |
|                19 | CH2CO           |                9 | CH2CO             | 1.7048 | 1.5542 |
|                20 | CHO             |               10 | CHO               | 0.7173 | 0.771  |
|                21 | CH3COO          |               11 | CCOO              | 1.27   | 1.6286 |
|                22 | CH2COO          |               11 | CCOO              | 1.27   | 1.4228 |
|                23 | HCOO            |               12 | HCOO              | 1.9    | 1.8    |
|                24 | CH3O            |               13 | CH2O              | 1.1434 | 1.6022 |
|                25 | CH2O            |               13 | CH2O              | 1.1434 | 1.2495 |
|                26 | CHO             |               13 | CH2O              | 1.1434 | 0.8968 |
|                27 | THF             |               43 | CY-CH2O           | 1.7023 | 1.8784 |
|                28 | CH3NH2          |               14 | CH2NH2            | 1.6607 | 1.6904 |
|                29 | CH2NH2          |               14 | CH2NH2            | 1.6607 | 1.3377 |
|                30 | CHNH2           |               14 | CH2NH2            | 1.6607 | 0.985  |
|                31 | CH3NH           |               15 | CH2NH             | 1.368  | 1.4332 |
|                32 | CH2NH           |               15 | CH2NH             | 1.368  | 1.0805 |
|                33 | CHNH            |               15 | CH2NH             | 1.368  | 0.7278 |
|                34 | CH3N            |               16 | (C)3N             | 1.0746 | 1.176  |
|                35 | CH2N            |               16 | (C)3N             | 1.0746 | 0.824  |
|                36 | ACNH2           |               17 | ACNH2             | 1.1849 | 0.8067 |
|                37 | AC2H2N          |               18 | PYRIDINE          | 1.4578 | 0.9022 |
|                38 | AC2HN           |               18 | PYRIDINE          | 1.2393 | 0.633  |
|                39 | AC2N            |               18 | PYRIDINE          | 1.0731 | 0.353  |
|                40 | CH3CN           |               19 | CH2CN             | 1.5575 | 1.5193 |
|                41 | CH2CN           |               19 | CH2CN             | 1.5575 | 1.1666 |
|                42 | COOH            |               20 | COOH              | 0.8    | 0.9215 |
|                43 | HCOOH           |               44 | HCOOH             | 0.8    | 1.2742 |
|                44 | CH2CL           |               21 | CCL               | 0.9919 | 1.3654 |
|                45 | CHCL            |               21 | CCL               | 0.9919 | 1.0127 |
|                46 | CCL             |               21 | CCL               | 0.9919 | 0.66   |
|                47 | CH2CL2          |               22 | CCL2              | 1.8    | 2.5    |
|                48 | CHCL2           |               22 | CCL2              | 1.8    | 2.1473 |
|                49 | CCL2            |               22 | CCL2              | 1.8    | 1.7946 |
|                50 | CHCL3           |               45 | CHCL3             | 2.45   | 2.8912 |
|                51 | CCL3            |               23 | CCL3              | 2.65   | 2.3778 |
|                52 | CCL4            |               24 | CCL4              | 2.618  | 3.1836 |
|                53 | ACCL            |               25 | ACCL              | 0.5365 | 0.3177 |
|                54 | CH3NO2          |               26 | CNO2              | 2.644  | 2.5    |
|                55 | CH2NO2          |               26 | CNO2              | 2.5    | 2.304  |
|                56 | CHNO2           |               26 | CNO2              | 2.887  | 2.241  |
|                57 | ACNO2           |               27 | ACNO2             | 0.4656 | 0.3589 |
|                58 | CS2             |               28 | CS2               | 1.24   | 1.068  |
|                59 | CH3SH           |               29 | CH3SH             | 1.289  | 1.762  |
|                60 | CH2SH           |               29 | CH3SH             | 1.535  | 1.316  |
|                61 | FURFURAL        |               30 | FURFURAL          | 1.299  | 1.289  |
|                62 | DOH             |               31 | DOH               | 2.088  | 2.4    |
|                63 | I               |               32 | I                 | 1.076  | 0.9169 |
|                64 | BR              |               33 | BR                | 1.209  | 1.4    |
|                65 | CH=-C           |               34 | C=-C              | 0.9214 | 1.3    |
|                66 | C=-C            |               34 | C=-C              | 1.303  | 1.132  |
|                67 | DMSO            |               35 | DMSO              | 3.6    | 2.692  |
|                68 | ACRY            |               36 | ACRY              | 1      | 0.92   |
|                69 | CL-(C=C)        |               37 | CLCC              | 0.5229 | 0.7391 |
|                70 | C=C             |                2 | C=C               | 1.2832 | 0.4582 |
|                71 | ACF             |               38 | ACF               | 0.8814 | 0.7269 |
|                72 | DMF             |               39 | DMF               | 2      | 2.093  |
|                73 | HCON(..         |               39 | DMF               | 2.381  | 1.522  |
|                74 | CF3             |               40 | CF2               | 1.284  | 1.266  |
|                75 | CF2             |               40 | CF2               | 1.284  | 1.098  |
|                76 | CF              |               40 | CF2               | 0.8215 | 0.5135 |
|                77 | COO             |               41 | COO               | 1.6    | 0.9    |
|                78 | CY-CH2          |               42 | CY-CH2            | 0.7136 | 0.8635 |
|                79 | CY-CH           |               42 | CY-CH2            | 0.3479 | 0.1071 |
|                80 | CY-C            |               42 | CY-CH2            | 0.347  | 0      |
|                81 | OH (S)          |                5 | OH                | 1.063  | 0.8663 |
|                82 | OH (T)          |                5 | OH                | 0.6895 | 0.8345 |
|                83 | CY-CH2O         |               43 | CY-CH2O           | 1.4046 | 1.4    |
|                84 | TRIOXAN         |               43 | CY-CH2O           | 1.0413 | 1.0116 |
|                85 | CNH2            |               14 | CH2NH2            | 1.6607 | 0.985  |
|                86 | NMP             |               46 | CY-CONC           | 3.981  | 3.2    |
|                87 | NEP             |               46 | CY-CONC           | 3.7543 | 2.892  |
|                88 | NIPP            |               46 | CY-CONC           | 3.5268 | 2.58   |
|                89 | NTBP            |               46 | CY-CONC           | 3.2994 | 2.352  |
|                91 | CONH2           |               47 | CONR              | 1.4515 | 1.248  |
|                92 | CONHCH3         |               47 | CONR              | 1.5    | 1.08   |
|                93 | HCONHCH3        |               49 | HCONR             | 2.4617 | 2.192  |
|                94 | HCONHCH2        |               49 | HCONR             | 2.4617 | 1.842  |
|               100 | CONHCH2         |               47 | CONR              | 1.5    | 1.08   |
|               101 | AM(CH3)2        |               48 | CONR2             | 2.4748 | 1.9643 |
|               102 | AMCH3CH2        |               48 | CONR2             | 2.2739 | 1.5754 |
|               103 | AM(CH2)2        |               48 | CONR2             | 2.0767 | 1.1866 |
|               104 | AC2H2S          |               52 | ACS               | 1.7943 | 1.34   |
|               105 | AC2HS           |               52 | ACS               | 1.6282 | 1.06   |
|               106 | AC2S            |               52 | ACS               | 1.4621 | 0.78   |
|               107 | H2COCH          |               53 | EPOXIDES          | 1.3601 | 1.8031 |
|               108 | COCH            |               53 | EPOXIDES          | 0.683  | 0.3418 |
|               109 | HCOCH           |               53 | EPOXIDES          | 0.9104 | 0.6538 |
|               110 | (CH2)2SU        |               56 | SULFONE           | 2.687  | 2.12   |
|               111 | CH2SUCH         |               56 | SULFONE           | 2.46   | 1.808  |
|               112 | (CH3)2CB        |               55 | CARBONAT          | 2.42   | 2.4976 |
|               113 | (CH2)2CB        |               55 | CARBONAT          | 2.42   | 2.0018 |
|               114 | CH2CH3CB        |               55 | CARBONAT          | 2.42   | 2.2497 |
|               119 | H2COCH2         |               53 | EPOXIDES          | 1.063  | 1.123  |
|               122 | CH3S            |               61 | CH2S              | 1.613  | 1.368  |
|               123 | CH2S            |               61 | CH2S              | 1.3863 | 1.06   |
|               124 | CHS             |               61 | CH2S              | 1.1589 | 0.748  |
|               153 | H2COC           |               53 | EPOXIDES          | 0.9104 | 0.6538 |
|               178 | C3H2N2+         |               84 | IMIDAZOL          | 1.3662 | 0.6797 |
|               179 | BTI-            |               85 | BTI               | 5.621  | 5.9463 |
|               184 | C3H3N2+         |               84 | IMIDAZOL          | 1.843  | 1.6997 |
|               189 | C4H8N+          |               87 | PYRROL            | 2.7867 | 2.7723 |
|               195 | BF4-            |               89 | BF4               | 3.9628 | 0.6214 |
|               196 | C5H5N+          |               90 | PYRIDIN           | 2.1094 | 2.5106 |
|               197 | OTF-            |               91 | OTF               | 3.371  | 2.0001 |
|               201 | -S-S-           |               93 | -S-S-             | 1.0678 | 2.244  |
|               209 | SO4             |               98 | SO4               | 0.9903 | 3.5249 |
|               210 | HSO4            |               98 | SO4               | 1.5654 | 3.8076 |
|               211 | PF6             |               99 | PF6               | 3.8183 | 3.6018 |
|               220 | C5H4N+          |               90 | PYRIDIN           | 2.4873 | 2.4457 |


## Using `ugropy` to retrieve UNIFAC-Dortmund subgroups

There is the possibility of using another library of our group
[`ugropy`](https://github.com/ipqa-research/ugropy) to retrieve the
UNIFAC-Dortmund subgroups and not suffer the pain of typing the subgroup
numbers and parameters by hand. The next Python snippet shows how you can use
it.

```python
from ugropy import dortmund, writers


names = ["water", "toluene", "cyclohexane"]
groups = [dortmund.get_groups(n).subgroups for n in names]

fortran_code = writers.to_yaeos(groups, dortmund)

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

molecules(3)%groups_ids = [78]
molecules(3)%number_of_groups = [6]
```

## Examples

Here is an example of a fully instantiated UNIFAC-Dortmund model. Please check
the `Gibbs Excess Models` section in the user documentation to learn all the
things you can do with this model.

Notice that here we are using the [[setup_dortmund]] function to instantiate
the model.

### Calculating activity coefficients
We can instantiate a [[UNIFAC]] model with a mixture ethanol-water and evaluate
the logarithm of activity coefficients of the model for a 0.5 mole fraction of
each, and a temperature of 298.0 K.

```fortran
use yaeos__constants, only: pr
use yaeos__models_ge_group_contribution_unifac, only: Groups, UNIFAC, setup_dortmund

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
model = setup_dortmund(molecules)

! Calculate ln_gammas
call model%ln_activity_coefficient([0.5_pr, 0.5_pr], 298.0_pr, ln_gammas)

print *, ln_gammas
```

## References

1. Weidlich, U., & Gmehling, J. (1987). A modified UNIFAC model. 1. Prediction
   of VLE, hE, and .gamma..infin. Industrial & Engineering Chemistry Research,
   26(7), 1372-1381. https://doi.org/10.1021/ie00067a018
2. Gmehling, J., Li, J., & Schiller, M. (1993). A modified UNIFAC model. 2.
   Present parameter matrix and results for different thermodynamic properties.
   Industrial & Engineering Chemistry Research, 32(1), 178-193.
   https://doi.org/10.1021/ie00013a024
3. Gmehling, J., Lohmann, J., Jakob, A., Li, J., & Joh, R. (1998). A Modified
   UNIFAC (Dortmund) Model. 3. Revision and Extension. Industrial & Engineering
   Chemistry Research, 37(12), 4876-4882. https://doi.org/10.1021/ie980347z
4. Lohmann, J., & Gmehling, J. (2001). Modified UNIFAC (Dortmund). Reliable
   Model for the Development of Thermal Separation Processes. JOURNAL OF
   CHEMICAL ENGINEERING OF JAPAN, 34(1), 43-54.
   https://doi.org/10.1252/jcej.34.43
5. Lohmann, J., Joh, R., & Gmehling, J. (2001). From UNIFAC to Modified UNIFAC
   (Dortmund). Industrial & Engineering Chemistry Research, 40(3), 957-964.
   https://doi.org/10.1021/ie0005710
6. Wittig, R., Lohmann, J., Joh, R., Horstmann, S., & Gmehling, J. (2001).
   Vapor−Liquid Equilibria and Enthalpies of Mixing in a Temperature Range from
   298.15 to 413.15 K for the Further Development of Modified UNIFAC
   (Dortmund). Industrial & Engineering Chemistry Research, 40(24), 5831-5838.
   https://doi.org/10.1021/ie010444j
7. Gmehling, J., Wittig, R., Lohmann, J., & Joh, R. (2002). A Modified UNIFAC
   (Dortmund) Model. 4. Revision and Extension. Industrial & Engineering
   Chemistry Research, 41(6), 1678-1688. https://doi.org/10.1021/ie0108043
8. Wittig, R., Lohmann, J., & Gmehling, J. (2003). Prediction of phase
   equilibria and excess properties for systems with sulfones. AIChE Journal,
   49(2), 530-537. https://doi.org/10.1002/aic.690490223
9. Jakob, A., Grensemann, H., Lohmann, J., & Gmehling, J. (2006). Further
   Development of Modified UNIFAC (Dortmund): Revision and Extension 5.
   Industrial & Engineering Chemistry Research, 45(23), 7924-7933.
   https://doi.org/10.1021/ie060355c
10. Hector, T., & Gmehling, J. (2014). Present status of the modified UNIFAC
    model for the prediction of phase equilibria and excess enthalpies for
    systems with ionic liquids. Fluid Phase Equilibria, 371, 82-92.
    https://doi.org/10.1016/j.fluid.2014.03.006
11. Constantinescu, D., & Gmehling, J. (2016). Further Development of Modified
    UNIFAC (Dortmund): Revision and Extension 6. Journal of Chemical &
    Engineering Data, 61(8), 2738-2748.
    https://doi.org/10.1021/acs.jced.6b00136
