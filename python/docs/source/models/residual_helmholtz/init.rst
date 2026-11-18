Residual Helmholtz Models
=========================

An equation of state relates temperature, volume, composition, and pressure.
Its residual Helmholtz energy measures departure from ideal-gas behavior;
derivatives provide properties and fugacity coefficients used in equilibrium
calculations.

Start with :doc:`cubic_eos/init` when you have critical constants and
acentric factors. :doc:`multifluid/init` provides a component-database model;
:doc:`saft/init` uses molecular segment parameters instead. These parameter
sets are not interchangeable.

All inherit the property and equilibrium methods in :doc:`/core`. Read the
:doc:`EoS tutorial </tutorial/equations_of_state/equations_of_state>` for a
calculation sequence, and distinguish a homogeneous volume root from a
stable multiphase state before interpreting density or energy results.

.. automodule:: yaeos.models.residual_helmholtz
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   cubic_eos/init
   multifluid/init
   saft/init
