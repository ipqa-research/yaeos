PC-SAFT
=======

``PCSAFT(m, sigma, epsilon_k, kij=None)`` describes each component by a
segment count ``m`` (dimensionless), segment diameter ``sigma`` (angstrom),
and energy scale ``epsilon_k`` (K). These are model parameters, not the
operating temperature or a molecular diameter to enter in litres.

Keep all three arrays and the square ``kij`` matrix in the composition's
component order. Omitting ``kij`` sets it to zero. ``set_kij(i, j, value)``
uses Python component indices starting at zero. This constructor does not
accept association-site or polar parameters; do not assume every published
PC-SAFT parameterization can be represented by this interface.

After construction, ``volume`` returns L for the supplied moles, not mass
density; divide total moles by volume for mol/L. For property units and
phase-result interpretation see :doc:`/core`. Continue with the
:doc:`SAFT tutorial </tutorial/equations_of_state/saft_eos/saft_eos>`.

.. automodule:: yaeos.models.residual_helmholtz.saft.pcsaft
    :members:
    :undoc-members:
    :show-inheritance:
