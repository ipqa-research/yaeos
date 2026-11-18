Multifluid Equations of State
=============================

Multifluid models combine component-specific descriptions with mixture
contributions. Here :doc:`gerg2008` selects components from a built-in set,
so construction uses names rather than arrays of critical constants.

Use the :doc:`multifluid tutorial </tutorial/equations_of_state/multifluid_eos/multifluid_eos>`
for an example and :doc:`/core` for inherited methods. A component being
available does not establish accuracy for every temperature, pressure,
composition, or phase: check the model's intended range for your application.

.. automodule:: yaeos.models.residual_helmholtz.multifluid
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   gerg2008
