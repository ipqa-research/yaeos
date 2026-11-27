.. include:: ../../README.md
   :parser: myst_parser.sphinx_

Start with a physical question
==============================

* **What pressure or volume does this fluid have?** Choose an
  :doc:`equation of state <models/residual_helmholtz/init>` and follow the
  :doc:`EoS tutorials <tutorial/equations_of_state/equations_of_state>`.
* **How nonideal is a liquid mixture?** Start with
  :doc:`excess Gibbs models <models/excess_gibbs/init>` for activity
  coefficients and liquid-liquid separation.
* **Will the mixture split into phases?** Read
  :doc:`phase equilibrium <tutorial/equations_of_state/phase_equilibrium/phase_equilibrium>`
  before using flash or saturation methods.

New users should follow the :doc:`tutorial/tutorial_tree` in order. The API
pages below give exact signatures; common property and equilibrium methods
are inherited from :doc:`core`, not repeated on every model page.

.. _python-input-conventions:

Input conventions to keep nearby
--------------------------------

The property methods use **kelvin**, **bar**, **litres**, and **moles**.
Energies reported in bar L convert to joules by multiplying by 100.
Check each signature: ``moles`` is an amount vector, whereas ``z`` in
equilibrium calculations is an overall mole-fraction vector summing to one.
A mole-fraction vector used as ``moles`` chooses a one-mole basis; the
corresponding volume is then the volume of that one mole.

Keep the same component order in every parameter array, composition vector,
and row/column of an interaction matrix. Outputs use that same order.
Mass fractions must be converted to mole fractions before use. Residual and
excess energies are contributions relative to a reference, not complete
enthalpies or Gibbs energies; see :doc:`core` before making an energy balance.


.. toctree::
   :maxdepth: 4
   :caption: Tutorial

   tutorial/tutorial_tree


.. toctree::
   :maxdepth: 2
   :caption: API documentation
   
   core

.. toctree::
   :maxdepth: 4

   models/init

.. toctree::
   :maxdepth: 5

   gpec
   envelopes
   fitting/init
   tools/init

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
