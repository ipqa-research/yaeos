Mixing Rules
============

A mixing rule combines pure-component EoS parameters into a mixture model.
It can change phase behavior substantially even when the pure-component
constants are unchanged. Construct the rule and pass it as ``mixrule`` to
the :doc:`cubic model <cubic_eos>`.

``QMR`` uses dimensionless binary ``kij`` and ``lij`` matrices. Each matrix
has one row and column per component in the same order as the critical
constants and composition vectors. ``QMRTD`` adds temperature dependence;
its reference temperatures are in K. Do not transfer fitted parameters to
a different mixing rule without checking their definition.

``HV`` and ``MHV`` couple an :doc:`excess Gibbs model </models/excess_gibbs/init>`
to the EoS. The liquid model and EoS must describe the same components in the
same order; MHV's ``q`` depends on the EoS. These rules do not mean that an
isolated excess-Gibbs model now supplies vapor properties on its own.

Follow the
:doc:`mixing-rule notebook </tutorial/equations_of_state/cubic_eos/cubic_mixing_rules>`
for construction examples, then :doc:`/fitting/init` if experimental data
are available to estimate interactions.

.. automodule:: yaeos.models.residual_helmholtz.cubic_eos.mixing_rules
    :members:
    :undoc-members:
    :show-inheritance:
