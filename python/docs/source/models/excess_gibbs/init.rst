Excess Gibbs Models
===================

These models describe how a liquid mixture differs from ideal mixing.
Their activity coefficients correct each component's liquid-phase chemical
potential: ``ln_gamma`` returns the natural logarithm, so use
``numpy.exp(model.ln_gamma(n, T))`` to obtain gamma. Gamma near one means a
small activity correction at that composition, not necessarily a stable
single liquid phase.

Use :doc:`nrtl` or :doc:`uniquac` with suitable fitted interaction parameters.
Use the UNIFAC variants when representing molecules by functional groups,
checking the variant's group definitions and available interactions first.
The variants are distinct parameterizations, not interchangeable names.

Inputs to property methods are component moles and temperature in K;
equilibrium methods take overall mole fractions. Keep the constructor's
component order throughout. ``flash_t`` describes a liquid-liquid split,
not boiling: these models alone do not provide a vapor equation of state
or pure-component vapor pressures. See :doc:`/core` for inherited methods
and :doc:`/tutorial/gibbs_excess/gibbs_excess` for examples.

.. automodule:: yaeos.models.excess_gibbs
   :undoc-members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1

   nrtl
   unifacvle
   unifac_dortmund
   psrk_unifac
   uniquac
