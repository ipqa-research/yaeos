UNIFAC Dortmund
===============

The Dortmund variant is a modified UNIFAC model with its own group and
interaction parameterization. Use ``UNIFACDortmund(molecules)`` when your
group decomposition and intended prediction match that variant; do not
substitute its parameters into the original UNIFAC model.

``molecules`` is a list of subgroup-count dictionaries, one per component
in the order later used by ``ln_gamma(n, T)``. The property call uses moles
and K and returns natural-log activity coefficients, not gamma directly.
Group counts are structural integers, not composition values.

Read :doc:`/models/groups`, the
:doc:`UNIFAC tutorial </tutorial/gibbs_excess/models/unifac_models>`, and
:doc:`/core` for excess-property units and liquid-liquid calculations.

.. automodule:: yaeos.models.excess_gibbs.dortmund
    :members:
    :undoc-members:
    :show-inheritance:
