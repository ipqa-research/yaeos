Equations of State (EoS)
========================

Use this route for pressure, density, and phase behavior of compressible
fluids. First construct a model from component parameters, then use **Bulk
properties** to evaluate a single homogeneous state. Only then move to
phase equilibrium, where phase compositions become additional unknowns.

The model chapters explain different parameter inputs: critical constants
for cubics, component names for GERG2008, and segment parameters for PC-SAFT.
They share the :doc:`/core` interface. Keep temperature in K, pressure in
bar, and volume in L; a property evaluated at a liquid root is not by itself
evidence that the feed remains entirely liquid.

After a single saturation point and flash are understood, try envelopes,
global binary diagrams, or fitting. Those tasks depend more strongly on
initial guesses and on the model's suitability for the mixture.


.. toctree::
   :maxdepth: 4

   bulk_properties.ipynb
   cubic_eos/cubic_eos.rst
   multifluid_eos/multifluid_eos.rst
   saft_eos/saft_eos.rst
   phase_equilibrium/phase_equilibrium.rst
   gpec/gpec.rst
   fitting/fitting.rst
