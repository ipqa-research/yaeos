EoS Fitting Module
==================

This package estimates binary equation-of-state interaction parameters from
experimental phase-equilibrium data. It is deliberately split into three
parts so that each piece can be understood on its own:

* :doc:`core` builds the model for a trial parameter vector and scores the
  model against the data.
* :doc:`model_setters` turn an optimizer's parameter vector into interaction
  parameters on the model.
* :doc:`solvers` compute the phase-equilibrium states that the objective
  function compares against the data.

Practical guidance
------------------

* Start with a model that already reproduces pure-component properties
  reasonably. Fitting interactions cannot repair a wrong pure-component
  description.
* Prepare the data as a pandas ``DataFrame`` with the documented columns
  (``kind``, ``x1``, ``y1``, ``T``, ``P``, optional ``weight``). Use
  kelvin, bar, and mole fractions, and keep the light component first.
* Not every data point needs every column. Leave unused values as ``NaN``;
  the fitter selects which columns matter from the ``kind`` label.
* The fitted parameters are specific to the model, the mixing rule, the
  chosen error definitions, and the data range. A good objective value does
  not guarantee a unique parameter set or valid extrapolation.

Work through the :doc:`fitting tutorial </tutorial/equations_of_state/fitting/fitting>`
for a complete example, then consult :doc:`core` for the exact fitter
signature.

.. automodule:: yaeos.fitting
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   core
   model_setters
   solvers