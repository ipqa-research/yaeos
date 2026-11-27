Alpha Functions
===============

The alpha function sets the temperature dependence of a cubic EoS
attractive parameter. Each cubic model has a default alpha function; an
alternative can be set on an existing model with ``model.set_alpha(...)``.

``AlphaSoave`` and ``AlphaRKPR`` take one dimensionless ``k`` value per
component, and ``AlphaMathiasCopeman`` takes three coefficient vectors
``c1``, ``c2`` and ``c3``. All vectors follow the same component order as the
critical constants. Changing the alpha function changes pure-component vapor
pressures, so check them before fitting binary interactions.

Follow the
:doc:`alpha-function notebook </tutorial/equations_of_state/cubic_eos/cubic_alphas>`
for examples.

.. automodule:: yaeos.models.residual_helmholtz.cubic_eos.alphas
    :members:
    :undoc-members:
    :show-inheritance:
