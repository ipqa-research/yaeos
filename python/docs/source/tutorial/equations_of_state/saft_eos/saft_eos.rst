SAFT Equations of State
=======================

This notebook introduces PC-SAFT construction and calculations. Prepare
segment numbers (dimensionless), diameters (angstrom), energy parameters
(K), and any binary interactions for the same parameterization. They are
not substitutes for the critical constants used by cubic EoS.

Keep component order fixed and distinguish parameter ``epsilon_k`` from
the operating temperature. Property calculations still use K, bar, L, and
mol. See :doc:`/models/residual_helmholtz/saft/pcsaft` for the constructor's
scope and :doc:`/core` for shared property and equilibrium methods.

.. toctree::
   :maxdepth: 4

   saft.ipynb
