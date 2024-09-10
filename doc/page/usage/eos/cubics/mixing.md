---
title: Mixing Rules
---

# Cubic Mixing Rules
All [[CubicEoS]] in `yaeos` include a [[CubicMixRule]] derived type, which 
handles how the \(D\) and \(B\) parameters in the CubicEoS are calculated.

## Quadratic Mixing Rules (QMR)
Quadratic mixing rules are the ussually most used mixing rules for cubic 
equations of state.

\[
    nB = \sum_i \sum_j n_i n_j \frac{b_i + b_j}{2} (1 - l_{ij})
\]

\[
    D = \sum_i \sum_j n_i n_j a_{ij}
\]

QMR are handled by the [[QMR]] derived type. Which can be used like:

```fortran
use yaeos, only: pr, QMR

type(QMR) :: mixrule
real(pr) :: kij(2, 2), lij(2, 2)

kij(1, :) = [0.0, 0.1]
kij(2, :) = [0.1, 0.0]

lij(1, :) = [0.0, 0.01]
lij(2, :) = [0.01, 0.0]

mixrule = QMR(k=kij, l=lij)
```

By default the \(a_{ij}\) matrix will be calculated with a constant \(k_{ij}\)
value (as shown below). But this functionality can be modified by replacing
the `get_aij` procedure

```fortran
use yaeos, only: pr, QMR

type(QMR) :: mixrule
real(pr) :: kij(2, 2), lij(2, 2)

kij(1, :) = [0.0, 0.1]
kij(2, :) = [0.1, 0.0]

lij(1, :) = [0.0, 0.01]
lij(2, :) = [0.01, 0.0]

mixrule = QMR(k=kij, l=lij)
mixrule%get_aij => my_aij_implementation

subroutine my_aij_implementation(self, ai, daidt, daidt2, aij, daijdt, daijdt2)
    class(QMR) :: self
    real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
    real(pr), intent(out) :: aij(:, :), daijdt(:, .), daijdt2(:, :)
    ! Implementation here
end subroutine
```

## \(G^E\) Models Mixing Rules
It is possible to mix the attractive parameter of Cubic Equations with an 
excess Gibbs-based model.
This can be useful for cases of polar molecules and/or systems that have been 
fitted to \(G^E\) models.

### Michelsen's Modified Huron-Vidal Mixing Rules
This mixing rule is based on the aproximate zero-pressure limit 
of a cubic equation of state. At the aproximate zero-pressure limit the
attractive parameter can be expressed as:

\[
\frac{D}{RTB}(n, T) = \sum_i n_i \frac{a_i(T)}{b_i} + \frac{1}{q}
\left(\frac{G^E(n, T)}{RT} + \sum_i n_i \ln \frac{B}{nb_i} \right)
\]

Where \(q\) is a weak function of temperature. In the case of `MHV`
and simplicity it is considered that depends on the model used.