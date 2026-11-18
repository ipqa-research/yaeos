Solvers Fitting module
======================

These helpers compute the phase-equilibrium states that the fitter compares
with experimental data. They answer a different question from the flash
methods in :doc:`core`: they return the **compositions of two phases at
given pressure and temperature**, not the phase amounts of an overall feed.

``solve_pt(model, pressure, temperature, kind)``
    Solves a binary equilibrium at the given pressure (bar) and temperature
    (K). It handles liquid-vapor and liquid-liquid cases; the ``kind``
    argument selects which (``"liquid-liquid"`` or ``"PT"``).
``binary_isofugacity_x1y1pt(x, p, t, model)``
    Evaluates the binary isofugacity condition. ``solve_pt`` does not call
    it (it runs ``flash_pt`` on a feed between the two initial guesses); it
    is available for building custom solvers. The argument
    ``x`` is the pair of trial mole fractions ``[y1, x1]`` of the light
    component, and the result is the squared residual of the equal-fugacity
    condition.
``find_init_binary_ll(model, pressure, temperature, kind)``
    Provides an initial guess for the liquid-liquid case by scanning the
    fugacity of the light component over the composition range.

Interpreting the output
-----------------------

The solver returns ``(x1, y1)``: the mole fraction of the **light component**
in the heavy and the light phase. The complete binary compositions are
``[x1, 1 - x1]`` and ``[y1, 1 - y1]`` in the model's component order. A
solution with ``x1 == y1`` means no separation was found at those
conditions, which the fitter interprets as a single-phase situation rather
than as an error.

Cautions
--------

Initial guesses and near-identical phase compositions make these solves
fragile. Check that the returned values are finite and physically plausible
(``0 <= x1, y1 <= 1``) and that the branch corresponds to the phase pair you
expect. A small isofugacity residual alone does not guarantee the correct
branch was followed.

These functions are used by :doc:`core` for the ``PT`` and ``liquid-liquid``
data kinds; see the
:doc:`fitting tutorial </tutorial/equations_of_state/fitting/fitting>` for
how they fit into a complete workflow.

.. automodule:: yaeos.fitting.solvers
    :members:
    :undoc-members:
    :show-inheritance:
