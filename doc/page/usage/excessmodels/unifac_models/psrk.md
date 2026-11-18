---
title: UNIFAC-PSRK
---

[TOC]

# Predictive Soave-Redlich-Kwong (PSRK) UNIFAC

## Model description

This page describes the **UNIFAC contribution used for PSRK**, not a complete
Predictive Soave–Redlich–Kwong equation of state. `setup_psrk(molecules)` returns
`type(UNIFAC)`, an excess-Gibbs liquid activity model. It does not construct an
SRK cubic, select its alpha function, attach a mixing rule, or calculate vapor
properties. Those choices must match the full PSRK formulation if that is your
goal; an arbitrary SRK/UNIFAC combination is not automatically PSRK.

This group set is relevant to predictive EOS work involving gases and solvents;
it includes dedicated light-gas groups absent from many original-UNIFAC examples.
For standalone liquid screening, compare the parameter coverage and validation
of [original UNIFAC](unifaclv.html) and [Dortmund](dortmund.html) rather than
selecting PSRK only because a subgroup appears in its table.

The model still combines geometric size/shape and energetic group-contact
contributions to *excess* Gibbs energy relative to an ideal liquid. Its internal
"residual" contribution is not an EOS ideal-gas departure. The returned
\(\gamma_i\) corrects liquid activity \(x_i\gamma_i\), not vapor pressure
directly and not the vapor fugacity coefficient. See the
[family overview](index.html) and [shared introduction](../index.html).

This is the Predictive Soave-Redlich-Kwong (PSRK) UNIFAC model. In this model,
the parameters are defined as:

$$
    z = 10
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

`T` is in K, \(a_{jk}\) is in K, \(b_{jk}\) is dimensionless, and \(c_{jk}\)
is in 1/K. The quadratic refers to \(U_{jk}/R\); the exponent is
\(-a_{jk}/T-b_{jk}-c_{jk}T\). Interaction coefficients are directional, so
reversing the indices can change the prediction. PSRK and Dortmund share this
functional form but **not** the same parameter set or combinatorial exponent.
The constructor loads the PSRK data; normally your inputs are just fragments
and their counts.

## Subgroups list

The list of the functional groups and its interaction parameters could be 
accessed on the DDBST web page:
[https://www.ddbst.com/psrk.html](https://www.ddbst.com/psrk.html)

We reproduce here the list of functional groups. To instantiate a UNIFAC-PSRK
model you must define which functional groups are used in a molecule by the
Subgroup Number column.

Pass subgroup IDs, not main-group IDs. `R` and `Q` are dimensionless relative
size and area. Dedicated whole-molecule groups such as CH4 and CO2 must not be
replaced by guessed atom-by-atom decompositions. The table is a published
reference; the installed library's bundled data determine actual availability.
Verify directional interaction coverage as well as subgroup existence.

These numbers belong to the **PSRK** variant only. For example, subgroup 118
is methane `CH4` here, while the same integer is the sulfone fragment
`(CH2)2SU` in the [original UNIFAC](unifaclv.html) table. Take subgroup
numbers from the table of the variant you are actually constructing.


|   Subgroup number | Subgroup name   | Main group    |      R |      Q |
|:-----------------:|:---------------:|:-------------:|:------:|:------:|
|                 1 | CH3             | [1] CH2       | 0.9011 | 0.848  |
|                 2 | CH2             | [1] CH2       | 0.6744 | 0.54   |
|                 3 | CH              | [1] CH2       | 0.4469 | 0.228  |
|                 4 | C               | [1] CH2       | 0.2195 | 0      |
|                 5 | CH2=CH          | [2] C=C       | 1.3454 | 1.176  |
|                 6 | CH=CH           | [2] C=C       | 1.1167 | 0.867  |
|                 7 | CH2=C           | [2] C=C       | 1.1173 | 0.988  |
|                 8 | CH=C            | [2] C=C       | 0.8886 | 0.676  |
|                 9 | ACH             | [3] ACH       | 0.5313 | 0.4    |
|                10 | AC              | [3] ACH       | 0.3652 | 0.12   |
|                11 | ACCH3           | [4] ACCH2     | 1.2663 | 0.968  |
|                12 | ACCH2           | [4] ACCH2     | 1.0396 | 0.66   |
|                13 | ACCH            | [4] ACCH2     | 0.8121 | 0.348  |
|                14 | OH              | [5] OH        | 1      | 1.2    |
|                15 | CH3OH           | [6] CH3OH     | 1.4311 | 1.432  |
|                16 | H2O             | [7] H2O       | 0.92   | 1.4    |
|                17 | ACOH            | [8] ACOH      | 0.8952 | 0.68   |
|                18 | CH3CO           | [9] CH2CO     | 1.6724 | 1.488  |
|                19 | CH2CO           | [9] CH2CO     | 1.4457 | 1.18   |
|                20 | HCO             | [10] HCO      | 0.998  | 0.948  |
|                21 | CH3COO          | [11] CCOO     | 1.9031 | 1.728  |
|                22 | CH2COO          | [11] CCOO     | 1.6764 | 1.42   |
|                23 | HCOO            | [12] HCOO     | 1.242  | 1.188  |
|                24 | CH3O            | [13] CH2O     | 1.145  | 1.088  |
|                25 | CH2O            | [13] CH2O     | 0.9183 | 0.78   |
|                26 | CHO             | [13] CH2O     | 0.6908 | 0.468  |
|                27 | THF             | [13] CH2O     | 0.9183 | 1.1    |
|                28 | CH3NH2          | [14] CNH2     | 1.5959 | 1.544  |
|                29 | CH2NH2          | [14] CNH2     | 1.3692 | 1.236  |
|                30 | CHNH2           | [14] CNH2     | 1.1417 | 0.924  |
|                31 | CH3NH           | [15] CNH      | 1.4337 | 1.244  |
|                32 | CH2NH           | [15] CNH      | 1.207  | 0.936  |
|                33 | CHNH            | [15] CNH      | 0.9795 | 0.624  |
|                34 | CH3N            | [16] (C)3N    | 1.1865 | 0.94   |
|                35 | CH2N            | [16] (C)3N    | 0.9597 | 0.632  |
|                36 | ACNH2           | [17] ACNH2    | 1.06   | 0.816  |
|                37 | C5H5N           | [18] PYRIDINE | 2.9993 | 2.113  |
|                38 | C5H4N           | [18] PYRIDINE | 2.8332 | 1.833  |
|                39 | C5H3N           | [18] PYRIDINE | 2.667  | 1.553  |
|                40 | CH3CN           | [19] CCN      | 1.8701 | 1.724  |
|                41 | CH2CN           | [19] CCN      | 1.6434 | 1.416  |
|                42 | COOH            | [20] COOH     | 1.3013 | 1.224  |
|                43 | HCOOH           | [20] COOH     | 1.528  | 1.532  |
|                44 | CH2CL           | [21] CCL      | 1.4654 | 1.264  |
|                45 | CHCL            | [21] CCL      | 1.238  | 0.952  |
|                46 | CCL             | [21] CCL      | 1.0106 | 0.724  |
|                47 | CH2CL2          | [22] CCL2     | 2.2564 | 1.988  |
|                48 | CHCL2           | [22] CCL2     | 2.0606 | 1.684  |
|                49 | CCL2            | [22] CCL2     | 1.8016 | 1.448  |
|                50 | CHCL3           | [23] CCL3     | 2.87   | 2.41   |
|                51 | CCL3            | [23] CCL3     | 2.6401 | 2.184  |
|                52 | CCL4            | [24] CCL4     | 3.39   | 2.91   |
|                53 | ACCL            | [25] ACCL     | 1.1562 | 0.844  |
|                54 | CH3NO2          | [26] CNO2     | 2.0086 | 1.868  |
|                55 | CH2NO2          | [26] CNO2     | 1.7818 | 1.56   |
|                56 | CHNO2           | [26] CNO2     | 1.5544 | 1.248  |
|                57 | ACNO2           | [27] ACNO2    | 1.4199 | 1.104  |
|                58 | CS2             | [28] CS2      | 2.057  | 1.65   |
|                59 | CH3SH           | [29] CH3SH    | 1.877  | 1.676  |
|                60 | CH2SH           | [29] CH3SH    | 1.651  | 1.368  |
|                61 | FURFURAL        | [30] FURFURAL | 3.168  | 2.484  |
|                62 | DOH             | [31] DOH      | 2.4088 | 2.248  |
|                63 | I               | [32] I        | 1.264  | 0.992  |
|                64 | BR              | [33] BR       | 0.9492 | 0.832  |
|                65 | CH=-C           | [34] C=-C     | 1.292  | 1.088  |
|                66 | C=-C            | [34] C=-C     | 1.0613 | 0.784  |
|                67 | DMSO            | [35] DMSO     | 2.8266 | 2.472  |
|                68 | ACRY            | [36] ACRY     | 2.3144 | 2.052  |
|                69 | CL-(C=C)        | [37] CLCC     | 0.791  | 0.724  |
|                70 | C=C             | [2] C=C       | 0.6605 | 0.485  |
|                71 | ACF             | [38] ACF      | 0.6948 | 0.524  |
|                72 | DMF             | [39] DMF      | 3.0856 | 2.736  |
|                73 | HCON(CH2)2      | [39] DMF      | 2.6322 | 2.12   |
|                74 | CF3             | [40] CF2      | 1.406  | 1.38   |
|                75 | CF2             | [40] CF2      | 1.0105 | 0.92   |
|                76 | CF              | [40] CF2      | 0.615  | 0.46   |
|                77 | COO             | [41] COO      | 1.38   | 1.2    |
|                78 | SIH3            | [42] SIH2     | 1.6035 | 1.2632 |
|                79 | SIH2            | [42] SIH2     | 1.4443 | 1.0063 |
|                80 | SIH             | [42] SIH2     | 1.2853 | 0.7494 |
|                81 | SI              | [42] SIH2     | 1.047  | 0.4099 |
|                82 | SIH2O           | [43] SIO      | 1.4838 | 1.0621 |
|                83 | SIHO            | [43] SIO      | 1.303  | 0.7639 |
|                84 | SIO             | [43] SIO      | 1.1044 | 0.4657 |
|                85 | NMP             | [44] NMP      | 3.981  | 3.2    |
|                86 | CCL3F           | [45] CCLF     | 3.0356 | 2.644  |
|                87 | CCL2F           | [45] CCLF     | 2.2287 | 1.916  |
|                88 | HCCL2F          | [45] CCLF     | 2.406  | 2.116  |
|                89 | HCCLF           | [45] CCLF     | 1.6493 | 1.416  |
|                90 | CCLF2           | [45] CCLF     | 1.8174 | 1.648  |
|                91 | HCCLF2          | [45] CCLF     | 1.967  | 1.828  |
|                92 | CCLF3           | [45] CCLF     | 2.1721 | 2.1    |
|                93 | CCL2F2          | [45] CCLF     | 2.6243 | 2.376  |
|                94 | AMH2            | [46] CON (AM) | 1.4515 | 1.248  |
|                95 | AMHCH3          | [46] CON (AM) | 2.1905 | 1.796  |
|                96 | AMHCH2          | [46] CON (AM) | 1.9637 | 1.488  |
|                97 | AM(CH3)2        | [46] CON (AM) | 2.8589 | 2.428  |
|                98 | AMCH3CH2        | [46] CON (AM) | 2.6322 | 2.12   |
|                99 | AM(CH2)2        | [46] CON (AM) | 2.4054 | 1.812  |
|               100 | C2H5O2          | [47] OCCOH    | 2.1226 | 1.904  |
|               101 | C2H4O2          | [47] OCCOH    | 1.8952 | 1.592  |
|               102 | CH3S            | [48] CH2S     | 1.613  | 1.368  |
|               103 | CH2S            | [48] CH2S     | 1.3863 | 1.06   |
|               104 | CHS             | [48] CH2S     | 1.1589 | 0.748  |
|               105 | MORPH           | [49] MORPH    | 3.474  | 2.796  |
|               106 | C4H4S           | [50] THIOPHEN | 2.8569 | 2.14   |
|               107 | C4H3S           | [50] THIOPHEN | 2.6908 | 1.86   |
|               108 | C4H2S           | [50] THIOPHEN | 2.5247 | 1.58   |
|               109 | H2C=CH2         | [2] C=C       | 1.3564 | 1.3098 |
|               110 | CH=-CH          | [34] C=-C     | 0.791  | 0.72   |
|               111 | NH3             | [55] NH3      | 0.851  | 0.778  |
|               112 | CO              | [63] CO       | 0.711  | 0.828  |
|               113 | H2              | [62] H2       | 0.416  | 0.571  |
|               114 | H2S             | [61] H2S      | 1.235  | 1.202  |
|               115 | N2              | [60] N2       | 0.856  | 0.93   |
|               116 | AR              | [59] AR       | 1.177  | 1.116  |
|               117 | CO2             | [56] CO2      | 1.3    | 0.982  |
|               118 | CH4             | [57] CH4      | 1.1292 | 1.124  |
|               119 | O2              | [58] O2       | 0.733  | 0.849  |
|               120 | D2              | [62] H2       | 0.37   | 0.527  |
|               121 | SO2             | [65] SO2      | 1.343  | 1.164  |
|               122 | NO              | [66] NO       | 0.716  | 0.62   |
|               123 | N2O             | [67] N2O      | 0.98   | 0.888  |
|               124 | SF6             | [68] SF6      | 2.374  | 2.056  |
|               125 | HE              | [69] HE       | 0.885  | 0.985  |
|               126 | NE              | [70] NE       | 0.886  | 0.986  |
|               127 | KR              | [71] KR       | 1.12   | 1.12   |
|               128 | XE              | [72] XE       | 1.13   | 1.13   |
|               129 | HF              | [73] HF       | 1.016  | 1.216  |
|               130 | HCL             | [74] HCL      | 1.056  | 1.256  |
|               131 | HBR             | [75] HBR      | 1.058  | 1.258  |
|               132 | HI              | [76] HI       | 1.393  | 1.208  |
|               133 | COS             | [77] COS      | 1.6785 | 1.316  |
|               134 | CHSH            | [29] CH3SH    | 1.425  | 1.06   |
|               135 | CSH             | [29] CH3SH    | 1.199  | 0.752  |
|               136 | H2COCH          | [51] EPOXY    | 1.3652 | 1.008  |
|               137 | HCOCH           | [51] EPOXY    | 1.1378 | 0.696  |
|               138 | HCOC            | [51] EPOXY    | 0.9104 | 0.468  |
|               139 | H2COCH2         | [51] EPOXY    | 1.5926 | 1.32   |
|               140 | H2COC           | [51] EPOXY    | 1.1378 | 0.78   |
|               141 | COC             | [51] EPOXY    | 0.6829 | 0.24   |
|               142 | F2              | [78] F2       | 0.75   | 0.88   |
|               143 | CL2             | [79] CL2      | 1.53   | 1.44   |
|               144 | BR2             | [80] BR2      | 1.9    | 1.66   |
|               145 | HCN             | [81] HCN      | 1.2    | 1.19   |
|               146 | NO2             | [82] NO2      | 1      | 1.1    |
|               147 | CF4             | [83] CF4      | 1.78   | 1.82   |
|               148 | O3              | [84] O3       | 1.1    | 1.27   |
|               149 | CLNO            | [85] CLNO     | 1.48   | 1.34   |
|               152 | CNH2            | [14] CNH2     | 0.9147 | 0.614  |


## Using `ugropy` to retrieve UNIFAC-PSRK subgroups

There is the possibility of using another library of our group
[`ugropy`](https://github.com/ipqa-research/ugropy) to retrieve the UNIFAC-PSRK
subgroups and not suffer the pain of typing the subgroup numbers and
parameters by hand. The next Python snippet shows how you can use it.

```python
from ugropy import psrk, writers


names = ["water", "toluene", "acetone"]
groups = [psrk.get_groups(n).subgroups for n in names]

fortran_code = writers.to_yaeos(groups, psrk)

print(fortran_code)
```

Typical output is a **setup snippet**, not a complete program. Confirm that
the assignments use PSRK subgroup definitions before copying them:

```fortran
use yaeos, only: Groups

type(Groups) :: molecules(3)

molecules(1)%groups_ids = [16]
molecules(1)%number_of_groups = [1]

molecules(2)%groups_ids = [9, 11]
molecules(2)%number_of_groups = [5, 1]

molecules(3)%groups_ids = [1, 18]
molecules(3)%number_of_groups = [1, 1]
```


## Examples

Here is an example of a fully instantiated UNIFAC-PSRK model. Please check the
`Gibbs Excess Models` section in the user documentation to learn all the things
you can do with this model.

Notice that here we are using the [[setup_psrk]] function to instantiate the
model.

### Calculating activity coefficients
We can instantiate a [[UNIFAC]] model with a mixture ethanol-water and evaluate
the logarithm of activity coefficients of the model for a 0.5 mole fraction of
each, and a temperature of 298.0 K.

```fortran
program psrk_unifac_activity
use yaeos, only: pr, Groups, UNIFAC, setup_psrk
implicit none

! Variables declarations
type(UNIFAC) :: model
type(Groups) :: molecules(2)
real(pr) :: ln_gammas(2)
real(pr) :: n(2)
real(pr) :: T
real(pr) :: Ge

! Variables instances
! Ethanol definition [CH3, CH2, OH]
molecules(1)%groups_ids = [1, 2, 14] ! Subgroups ids
molecules(1)%number_of_groups = [1, 1, 1] ! Subgroups occurrences

! Water definition [H2O]
molecules(2)%groups_ids = [16]
molecules(2)%number_of_groups = [1]

! Model setup
model = setup_psrk(molecules)

! Calculate ln_gammas
n = [0.5_pr, 0.5_pr]
T = 298.0_pr
call model%ln_activity_coefficient(n, T, ln_gammas)
call model%excess_gibbs(n, T, Ge=Ge)

print *, "ln(gamma): ethanol, water = ", ln_gammas
print *, "gamma = ", exp(ln_gammas)
print *, "Molar GE [bar L/mol] = ", Ge/sum(n)
end program psrk_unifac_activity
```

This complete program deliberately evaluates only the group model. Output order
is ethanol then water, matching `molecules` and `n`. Counts describe fragments
per molecule, not mixture composition. `ln_gammas` are natural logarithms; their
exponentials are dimensionless activity coefficients. The molar excess energy
is in bar L/mol (multiply by 100 for J/mol), not the full mixing Gibbs energy.
No vapor phase or gas solubility is solved here.

Do not assume this program prints the same activity coefficients as original
UNIFAC. Even if some parameter entries coincide, geometry, group definitions,
interaction tables, and temperature handling are variant-specific. Compare
numerical results only after checking the installed parameter tables and the
exact groups and conventions used by both constructors.

## Checks, common errors, and next steps

Use equally sized subgroup-ID/count arrays, list each subgroup once with a
positive total count, and check that all required directed interactions are
available. Do not use an apparent zero or missing table entry as proof of ideal
behavior. Keep the PSRK data, temperature law, and geometry together; borrowed
coefficients from other variants invalidate the model definition. In
particular, subgroup numbers must come from the PSRK table on this page, since
the same integers mean different fragments in the other variants.

Use kelvin and nonnegative mole amounts with positive total. Scaling all amounts
should not change `gamma`; the present pure component should have a coefficient
approaching 1. These checks test implementation and input consistency, not
accuracy.

Be clear about what the PSRK name implies and what it does not. High-pressure
reliability of a full PSRK equation of state does not follow from evaluating
this pressure-independent group model, especially for supercritical light gases
where an ordinary pure-liquid standard state needs careful treatment. If your
goal is a high-pressure vapor-liquid calculation, you need the equation of
state itself.

For thermal properties see the [shared GeModel methods](../index.html). For EOS
coupling, continue to the [cubic mixing rules](../../eos/cubics/mixing.html) and
use a consistent full formulation. Always validate phase-equilibrium predictions
against data in the intended range.

## References

1. Holderbaum T., "Die Vorausberechnung von Dampf-Flüssig-Gleichgewichten mit
   einer Gruppenbeitragszustandsgleichung", Thesis, Universität Dortmund, 1990
2. Holderbaum T., Gmehling J., "PSRK: Eine Zustandsgleichung zur Vorhersage von
   Dampf/Flüssig- Gleichgewichten bei mittleren und hohen Drücken.",
   Chem.Ing.Tech. CIT, 63(1), 57-59, 1991
3. Holderbaum T., Gmehling J., "PSRK: A Group-Contribution Equation of State
   based on UNIFAC", Fluid Phase Equilib., 70, 251-265, 1991
4. Fischer K., "Die PSRK-Methode: Eine Zustandgleichung unter Verwendung des
   UNIFAC-Gruppenbeitragsmodells", Thesis, C.-v.-O. Universität Oldenburg, 1993
5. Fischer K., Gmehling J., "Further Development, Status and Results of the
   PSRK Method for the Prediction of Vapor-Liquid Equilibria and Gas
   Solubilities", Fluid Phase Equilib., 112, 1-22, 1995
6. Fischer K., Gmehling J., "Further Development, Status and Results of the
   PSRK Method for the Prediction of Vapor-Liquid Equilibria and Gas
   Solubilities", Fluid Phase Equilib., 121, 185-206, 1996
7. Gmehling J., Li J., Fischer K., "Further development of the PSRK model for
   the prediction of gas solubilities and vapor-liquid-equilibria at low and
   high pressures", Fluid Phase Equilib., 141, 113-127, 1997
8. Horstmann S., Fischer K., Gmehling J., "PSRK group contribution equation of
   state: revision and extension III", Fluid Phase Equilib., 167, 173-186, 2000
9. Horstmann S., "Theoretische und experimentelle Untersuchungen zum
   Hochdruckphasengleichgewichtsverhalten fluider Stoffgemische für die
   Erweiterung der PSRK-Gruppenbeitragszustandsgleichung", Thesis, C.-v.-O.
   Universität Oldenburg, 2000
10. Horstmann S., Jabloniec A., Krafczyk J., Fischer K., Gmehling J., "PSRK
    Group Contribution Equation of State: Comprehensive Revision and Extension
    IV, Including Critical Constants and a-Function Parameters for 1000
    Components", Fluid Phase Equilib., 227(2), 157-164, 2005
