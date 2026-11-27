UNIQUAC
=======

UNIQUAC combines molecular size/shape effects with energetic interactions
to describe liquid nonideality. Supply relative surface areas ``qs`` and
relative volumes ``rs`` as dimensionless arrays in component order. Notice
that the constructor takes **q before r**.

Interaction matrices use the same component order on both axes and are
directional, not generally symmetric. The implemented temperature dependence
is

.. math::

   \tau_{lk} = \exp\left(a_{lk} + \frac{b_{lk}}{T}
   + c_{lk}\ln T + d_{lk} T + e_{lk} T^2\right)

so ``b`` carries the usual energy scale expressed in kelvin (it is
:math:`-\Delta U_{lk}/R`), and the other terms are optional corrections that
are zero unless you supply them. Check the temperature dependence, sign, and
energy-to-temperature conversion of any published parameter set before
assigning it to ``aij`` through ``eij``. Omitting an interaction matrix
supplies zeros, not estimated interactions.

Evaluate ``ln_gamma(n, T)`` with component moles and temperature in K; its
exponential gives dimensionless activity coefficients. Continue with the
:doc:`UNIQUAC notebook </tutorial/gibbs_excess/models/uniquac>` and
:doc:`/core` for liquid-liquid equilibrium and excess properties.

.. automodule:: yaeos.models.excess_gibbs.uniquac
    :members:
    :undoc-members:
    :show-inheritance:
