Core Fitting module
===================

``BinaryFitter`` scores a model against experimental data and drives an
optimizer over the interaction parameters. It needs three things from you:
a **model setter** that builds a model from a parameter vector, the **data**
to fit, and (optionally) the error definitions to use.

Data format
-----------

The data is a pandas ``DataFrame``. The documented columns are:

``kind``
    Which measurement was made. The documented values are ``bubbleT``,
    ``bubbleP``, ``dewT``, ``dewP``, ``PT``, ``liquid-liquid``, and
    ``critical``; the shorter aliases ``bubble`` and ``dew`` are also
    accepted. Use ``critical`` rather than the ``CP`` alias: the critical
    line is only computed when at least one row has ``kind == "critical"``,
    so a data set that marks critical points only as ``CP`` fails.
``T``
    Temperature in K.
``P``
    Pressure in bar.
``x1``
    Mole fraction of component 1 (the light component) in the heavy phase.
``y1``
    Mole fraction of component 1 (the light component) in the light phase.
``weight`` (optional)
    Relative weight of the point; points without a weight count as 1.

What each ``kind`` does with the data:

* ``bubbleP`` compares the experimental pressure with the bubble saturation
  pressure at ``T`` for the liquid of composition ``x1``.
* ``bubbleT`` compares the experimental temperature with the bubble saturation
  temperature at ``P`` for the liquid of composition ``x1``.
* ``dewP`` and ``dewT`` do the same for the dew point, using the vapor
  composition ``y1``.
* ``PT`` and ``liquid-liquid`` compare the two-phase equilibrium computed at
  ``P`` and ``T`` with the experimental compositions.
* ``critical`` compares the experimental point with the nearest point on the
  model's calculated critical line.

``x1`` and ``y1`` are mole fractions of the *same* component in the two
phases, not fractions of two different components. Columns that a given
``kind`` does not use may be left as ``NaN``. An unrecognized ``kind`` raises
``ValueError``.

How the fit runs
----------------

``fit(x0, bounds=None, method="Nelder-Mead", optimizer_options=None)``
minimizes the objective with ``scipy.optimize.minimize``. ``x0`` is the
initial parameter guess, ``bounds`` may be ``None`` or a sequence of
``(min, max)`` pairs, one per fitted parameter, and ``optimizer_options`` is
passed to ``minimize`` as ``options``. The optimizer result is stored in the
``solution`` property, and the ``model`` property rebuilds the fitted model
as ``model_setter(solution.x, *model_setter_args)``.

Errors and diagnostics
----------------------

The default pressure error is the squared difference between the
experimental and model pressures divided by the experimental pressure; the
temperature error does the same with temperatures. Points at low values
therefore weigh more heavily than the raw difference suggests. The default
composition error is the sum of the absolute log ratios of the two
light-component mole fractions, which punishes near-zero compositions
heavily. Supplying your own ``pressure_error``, ``temperature_error``, or
``composition_error`` callables changes what the objective value means, so
keep the definitions consistent when comparing fits.

Read the residuals stored in ``evaluations["residuals"]`` and the per-point
contributions in ``evaluations["contributions"]``: an overall objective value
can hide a few badly predicted points. Check the fitted model against data
that were not used in the fit, and confirm that the phases found are the
expected ones before trusting the parameters.

See :doc:`model_setters` for the parameter-vector conventions and
:doc:`solvers` for the equilibrium calculations used by the ``PT`` and
``liquid-liquid`` data kinds.

.. automodule:: yaeos.fitting.core
    :members:
    :undoc-members:
    :show-inheritance:
