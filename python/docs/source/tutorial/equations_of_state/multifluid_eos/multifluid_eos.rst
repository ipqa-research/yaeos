Multifluid Equations of State
=============================

Use the GERG2008 notebook to build a mixture from supported component names
and evaluate it through the common EoS interface. Unlike cubic models,
this constructor selects built-in component parameters rather than asking
you for critical constants.

The names' order defines the order of every composition vector. From a
returned volume in L, ``sum(n) / V`` gives molar density in mol/L, not mass
density. Read :doc:`/models/residual_helmholtz/multifluid/gerg2008` and
:doc:`/core` before treating residual energies or heat capacities as total
properties.

.. toctree::
   :maxdepth: 4

   gerg2008.ipynb
