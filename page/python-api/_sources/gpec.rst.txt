GPEC module
===========

``GPEC`` is a driver for global phase-diagram calculations of **binary**
mixtures. Given an existing equation-of-state model, it traces critical lines
and three-phase lines, and provides isothermal (Pxy) and isobaric (Txy)
views of the mixture's phase behavior. This is a broader task than
a flash at one feed composition.

Before constructing it
----------------------

* Build and check a two-component model first: verify a property, a
  saturation point, and a flash. GPEC inherits every weakness of that model.
* Keep the component order fixed. The driver works with the first two
  components of the model, and its "first" and "second" component references
  (almost pure starting compositions, critical-line directions) follow that
  order.
* ``max_pressure`` is in bar and bounds the tracing range; ``max_points``
  limits the number of traced points, so a truncated curve usually means a
  numerical limit rather than a physical endpoint.

Reading the results
-------------------

``plot_gped`` draws the pure saturation curves, the critical line, and any
three-phase lines it found. ``calc_pxy`` and ``calc_txy`` return the
underlying phase-equilibrium data used for those views; ``plot_pxy`` and
``plot_txy`` add the plot. Composition coordinates are mole fractions in the
model's fixed order.

Continuation steps (``step_12``, ``step_21``), starting compositions
(``x10``, ``x20``), and the ``stability_analysis`` switch influence which
branches are discovered. A missing branch is not proof that it does not
exist, and a traced branch is not proof that it is globally stable: check
against stability analysis and, where relevant, experimental data.

Read the :doc:`GPEC tutorial <tutorial/equations_of_state/gpec/gpec>` for a
worked example, :doc:`core` for the underlying equilibrium methods, and
:doc:`envelopes` for the data structures used by envelope calculations.

.. automodule:: yaeos.gpec
    :members:
    :undoc-members:
    :show-inheritance:
