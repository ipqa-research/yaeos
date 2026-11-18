NRTL
====

NRTL describes nonideal liquid mixing through fitted binary interactions.
Use it for activity coefficients or liquid-liquid separation when parameters
are available for the components and temperature range of interest.

The implemented convention is ``tau_ij = a_ij + b_ij / T`` and
``G_ij = exp(-c_ij * tau_ij)``. Thus ``a`` and ``c`` are dimensionless and
``b`` is in K. Each is a square component-by-component matrix. Directional
interactions need not be symmetric: do not replace ``a12`` by ``a21`` or
transpose a literature table without checking its convention. Energy-valued
NRTL parameters require conversion to this convention before use.

``ln_gamma(n, T)`` returns one natural-log activity coefficient per component
in the constructor's order; use ``numpy.exp`` for gamma. See
:doc:`/core` for derivatives and liquid-liquid flashes, and the
:doc:`NRTL notebook </tutorial/gibbs_excess/models/nrtl>` for a worked setup.

.. automodule:: yaeos.models.excess_gibbs.nrtl
    :members:
    :undoc-members:
    :show-inheritance:
