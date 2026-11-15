Cubic Equations of State
========================

Read these notebooks in order: choose a cubic model, examine its
temperature-dependent attraction (the alpha function), then specify how
components mix. You need critical temperatures in K, critical pressures in
bar, and dimensionless acentric factors in a consistent component order.

Changing an alpha function changes pure-fluid temperature dependence;
changing a mixing rule changes mixture interactions. Do not use binary
interaction parameters to hide an incorrect unit conversion or component
ordering. The :doc:`API model pages </models/residual_helmholtz/cubic_eos/init>`
give constructor signatures; :doc:`/core` provides inherited calculations.

.. toctree::
   :maxdepth: 4

   cubic_models.ipynb
   cubic_alphas.ipynb
   cubic_mixing_rules.ipynb
