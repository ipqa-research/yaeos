GPEC (Global Phase Equilibrium Calculations)
=============================================

GPEC explores the organization of a **binary** mixture's phase diagram,
including critical lines and connections to three-phase boundaries. Use it
when a single flash or fixed-composition envelope does not explain the
mixture's overall phase behavior.

First complete the
:doc:`phase-equilibrium tutorial <../phase_equilibrium/phase_equilibrium>`
with a two-component EoS. Check that its parameters are appropriate over
the wider conditions being traced. Pressure limits use bar; composition
coordinates are mole fractions in the model's fixed component order.

Trace limits and continuation step sizes affect which branches are found;
a missing branch is not proof that it cannot exist. Inspect the calculated
critical and three-phase information, not just the plotted lines. The
:doc:`/gpec` reference describes the driver and its controls.

.. toctree::
   :maxdepth: 4

   gpec_tutorial.ipynb
