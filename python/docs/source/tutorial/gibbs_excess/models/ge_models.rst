Models
======

Choose NRTL or UNIQUAC when you have fitted interaction parameters in the
implemented convention. Choose a UNIFAC variant when you have a compatible
functional-group decomposition and the needed interactions are available.
Prediction from groups still depends on the selected parameterization and
its range of applicability.

The notebooks below explain construction. Before copying literature
parameters, check matrix direction, temperature dependence, units, and
component order. NRTL and UNIQUAC interaction matrices generally need not
be symmetric. UNIFAC dictionaries contain subgroup counts, not mole
fractions. Return to :doc:`/models/excess_gibbs/init` for API links and
:doc:`/core` for activity coefficients, derivatives, and liquid splits.


.. toctree::
   :maxdepth: 4

   nrtl.ipynb
   uniquac.ipynb
   unifac_models.ipynb
