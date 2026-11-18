SAFT Equations of State
=======================

SAFT models describe fluids using molecular segments rather than only
critical constants. The Python model here is :doc:`pcsaft`; use it when you
have a compatible set of segment parameters for your components.

It inherits the same :doc:`/core` property and equilibrium methods as cubic
EoS. Begin with the :doc:`SAFT tutorial </tutorial/equations_of_state/saft_eos/saft_eos>`
and verify a homogeneous property before attempting saturation or a flash.

.. automodule:: yaeos.models.residual_helmholtz.saft
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   pcsaft
