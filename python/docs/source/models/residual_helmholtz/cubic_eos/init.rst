Cubic Equations of State
========================

Cubic EoS are a practical starting point for compressible-fluid properties
and vapor-liquid equilibrium. Supply one critical temperature (K), critical
pressure (bar), and dimensionless acentric factor per component. The mixing
rule then describes how unlike components interact.

Start with :doc:`cubic_eos`, then read :doc:`mixing_rules` before adjusting
binary interactions. Pure-component constants alone do not guarantee good
mixture predictions. The
:doc:`cubic tutorial </tutorial/equations_of_state/cubic_eos/cubic_eos>`
also introduces :doc:`alpha functions <alphas>`, which control temperature
dependence of attraction. Shared calculations are documented in :doc:`/core`.

.. automodule:: yaeos.models.residual_helmholtz.cubic_eos
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   cubic_eos
   alphas
   mixing_rules
