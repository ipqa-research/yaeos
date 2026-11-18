Model Setters Fitting module
============================

A model setter converts the optimizer's trial parameter vector into a model
that the :doc:`core` fitter can score. The convention used throughout the
package is ``setter(parameters, *args) -> ArModel``: the parameters come
first, and the extra arguments are the ones passed to the fitter as
``model_setter_args``.

The supplied setters
--------------------

``fit_kij_lij(x, model, fit_kij, fit_lij)``
    Updates the model's mixing rule with binary interaction parameters. When
    both are fitted, the vector is ``[kij, lij]``; when only one is fitted,
    the vector has a single element. The matrices are built symmetric, so
    this setter cannot represent directional 12/21 differences.
``fit_mhv_nrtl(x, args)`` and ``fit_hv_nrtl(x, args)``
    Build an NRTL excess-Gibbs model and couple it to the cubic EoS through
    the MHV or HV mixing rule. Their extra arguments arrive as a single
    tuple, ``args = (model, q)`` for MHV and ``args = (model,)`` for HV, so
    pass them to the fitter wrapped once more, for example
    ``model_setter_args=((model, q),)``. The parameter vector is
    ``[a12, a21, b12, b21, alpha]``, matching the NRTL convention
    ``tau_ij = a_ij + b_ij / T`` with ``G_ij = exp(-c_ij tau_ij)``: ``a`` and
    ``alpha`` are dimensionless, ``b`` is in K, and the 12 and 21 entries are
    independent.

Practical notes
---------------

* The vector order is part of the interface. A vector of the right length
  but the wrong order will silently fit the wrong parameters.
* Setters mutate the model and return it. Keep the parameter vector yourself
  if you need to know which values produced a given result.
* Writing your own setter is normal and often better than using the
  supplied ones: it lets you fit parameters the helpers do not cover, such
  as asymmetric ``kij`` values or several parameters at once. Follow the
  same calling convention and return the updated model.

See :doc:`core` for how the setter is called and the
:doc:`fitting tutorial </tutorial/equations_of_state/fitting/fitting>` for a
complete worked example, including a custom setter.

.. automodule:: yaeos.fitting.model_setters
    :members:
    :undoc-members:
    :show-inheritance:
