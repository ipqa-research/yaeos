UNIFAC PSRK
===========

``UNIFACPSRK`` supplies the excess-Gibbs group-contribution model used in
the PSRK approach. It is not itself the ``PSRK`` equation of state; use
the :doc:`cubic EoS constructor </models/residual_helmholtz/cubic_eos/cubic_eos>`
when pressure, volume, or vapor properties are required.

Supply one subgroup-count dictionary per molecule using the PSRK group's
definitions and available interactions. Keep the list order for subsequent
mole vectors. ``ln_gamma(n, T)`` takes mol and K and returns natural logs
of dimensionless activity coefficients.

See :doc:`/models/groups` and
:doc:`/tutorial/gibbs_excess/models/unifac_models` for the input workflow;
:doc:`/core` documents the inherited excess-Gibbs methods.

.. automodule:: yaeos.models.excess_gibbs.psrk_unifac
    :members:
    :undoc-members:
    :show-inheritance:
