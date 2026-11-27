UNIFACVLE
=========

This is the original UNIFAC VLE parameterization for estimating liquid
activity coefficients from functional groups. Construct ``UNIFACVLE`` with
one subgroup-count dictionary per molecule; see :doc:`/models/groups`.
Counts describe a molecule's structure, not the feed composition.

For example, ``UNIFACVLE([{16: 1}, {1: 1, 2: 1, 14: 1}])`` describes
water followed by ethanol. A later call
``model.ln_gamma([0.8, 0.2], 298.15)`` uses that same component order,
takes mole amounts (here one mole in total) and a temperature in kelvin, and
returns the natural logarithms of the activity coefficients.

Check available subgroup interactions and the parameterization's intended
use; the name VLE does not mean the object alone calculates vapor pressure.
Continue with :doc:`/tutorial/gibbs_excess/models/unifac_models` and
:doc:`/core` for the shared liquid-property interface.

.. automodule:: yaeos.models.excess_gibbs.unifac
    :members:
    :undoc-members:
    :show-inheritance:
