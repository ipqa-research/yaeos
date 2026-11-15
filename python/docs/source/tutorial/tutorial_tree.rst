Tutorial
========

Suggested route
---------------

1. Open **Getting started** to learn how the Python interface relates to the
   Fortran library and where inherited methods are documented.
2. Choose **Equations of state** for pressure, density, vapor-liquid
   equilibrium, and phase envelopes. Begin with a cubic model and a bulk
   property calculation before attempting a phase split.
3. Choose **Excess Gibbs models** for liquid activity coefficients and
   liquid-liquid equilibrium. These models alone do not supply vapor
   pressures or gas-phase densities.
4. Attempt fitting or global phase diagrams only after reproducing a simple
   property and equilibrium calculation with a fixed model.

The notebooks below contain worked calculations. Keep
:ref:`python-input-conventions` open while adapting them: changing the
component order or using Celsius instead of kelvin changes the physical
problem, even when the code still runs. For exact method arguments and
return values, consult :doc:`/core` alongside the notebook.

.. toctree::
   :maxdepth: 4

   getting_started.ipynb
   equations_of_state/equations_of_state.rst
   gibbs_excess/gibbs_excess.rst
