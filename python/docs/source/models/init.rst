Models
======

Choose the model before the solver
----------------------------------

A model describes the fluid's nonideal behavior; an equilibrium solver asks
whether that description predicts one phase or several at specified
conditions. Solver convergence cannot compensate for unsuitable parameters.

* :doc:`Residual Helmholtz models <residual_helmholtz/init>` are equations of
  state: use them when pressure, volume, or vapor-liquid equilibrium matters.
* :doc:`Excess Gibbs models <excess_gibbs/init>` describe departures from
  ideal liquid mixing: use them for activity coefficients and liquid-liquid
  separation, or as part of an EoS mixing rule.
* :doc:`groups` explains the group dictionaries used by UNIFAC models.

Choose parameters for the actual components and temperature range, keep
their order fixed, and compare predictions with data relevant to your task.
The common property and phase-equilibrium methods are documented in
:doc:`/core`; the pages below document model construction.

.. automodule:: yaeos.models
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   residual_helmholtz/init
   excess_gibbs/init
   groups
