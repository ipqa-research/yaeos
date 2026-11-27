Excess Gibbs models
===================

Use this route for nonideal liquids. Start with **Bulk properties** to
understand activity coefficients and excess energies, choose a model in
**Models**, then work through **Equilibrium calculations** for liquid-liquid
separation.

The central quantity ``ln_gamma(n, T)`` is the natural logarithm of the
activity coefficient. Its exponential is dimensionless; it is not a mole
fraction. Temperature is in K and the amount vector is in mol, always in
the model's component order. Excess energies are in bar L for that amount
basis, not automatically J/mol.

These models describe liquid mixing, not a gas equation of state. A
two-phase ``flash_t`` result labels two liquids, not liquid and vapor.
For vapor-liquid calculations involving pressure and density, continue
with :doc:`/tutorial/equations_of_state/equations_of_state`.
See :doc:`/core` for inherited method signatures.


.. toctree::
   :maxdepth: 4

   bulk_properties.ipynb
   models/ge_models.rst
   equilibrium_calcs.ipynb

