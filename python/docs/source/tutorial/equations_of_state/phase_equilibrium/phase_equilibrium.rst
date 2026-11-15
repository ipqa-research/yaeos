Phase Equilibrium Calculations
==============================

Which question are you asking?
------------------------------

* **At what pressure does boiling begin?** Use
  ``saturation_pressure(z, temperature, kind="bubble")``. The specified
  composition is the liquid; the returned vapor is just appearing.
* **At what pressure does condensation begin?** Use ``kind="dew"``. The
  specified composition is the vapor; the returned liquid is just appearing.
* **How much of each phase exists at a given pressure and temperature?** Use
  ``flash_pt(z, pressure, temperature)`` with the overall feed composition.
* **Where does a phase boundary go as conditions change?** Trace an envelope
  after understanding and checking a single saturation point.

Prerequisites and interpretation
--------------------------------

Construct an :doc:`EoS model </models/residual_helmholtz/init>` first and
verify a bulk-property calculation. Supply normalized mole fractions in the
model's component order, temperature in K, and pressure in bar. Saturation
pressure guesses ``p0`` are also in bar; ``y0`` is an initial composition
guess for the incipient phase, not an imposed second-phase composition.

Two-phase results contain ``x`` (heavy phase), ``y`` (light phase), and
``beta`` (light-phase mole fraction). A saturation result has the incipient
phase with zero amount: ``beta`` is 0 at a bubble point and 1 at a dew point.
Check phase normalization and
``z = (1 - beta) * x + beta * y`` before using a flash result. The documented
``flash_pt`` failure value is ``beta == -1``; a failed saturation calculation
reports a non-finite pressure instead. A single-phase state is not an error
merely because two distinct phase compositions cannot be found.

Near critical points the phases become similar, making initialization and
interpretation harder. Stability analysis helps distinguish physical
separation from an unhelpful or metastable solution; convergence alone is
not a stability test. See :doc:`/core` for exact arguments and
:doc:`/envelopes` for the arrays returned by envelope calculations.

.. toctree::
   :maxdepth: 4

   phase_equilibrium_tutorial.ipynb
