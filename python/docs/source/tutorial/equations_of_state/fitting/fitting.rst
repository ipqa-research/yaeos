Fitting Equations of State
==========================

Use fitting to estimate binary interaction parameters from experimental
phase-equilibrium data, after you can calculate those equilibria with a
fixed model. The workflow is: prepare consistent data, choose which
parameters may vary, minimize an error measure, and check predictions
against data not used for fitting.

Use K, bar, and mole fractions, and keep the light component first as
required by ``BinaryFitter``. Check phase labels and experimental
uncertainties before combining datasets. A small objective value does not
establish unique parameters or reliable extrapolation.

Read :doc:`/fitting/init` for the fitter, model setters, and equilibrium
solver roles. Bounds, initial guesses, and the chosen error functions are
part of the modeling decision, not just optimizer settings.

.. toctree::
   :maxdepth: 4

   fitting_tutorial.ipynb
