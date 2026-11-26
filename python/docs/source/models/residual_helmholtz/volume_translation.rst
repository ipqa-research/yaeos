Volume translations
===================

Constant translation
--------------------

``ConstantVolumeTranslation`` wraps an existing residual Helmholtz model.
Provide one constant molar shift per component, in the same order as the base
model. The model copies the configured base EOS when constructed.

.. code-block:: python

   from yaeos import ConstantVolumeTranslation, PengRobinson76

   base = PengRobinson76([190.564, 305.32], [45.99, 48.72], [0.01142, 0.0995])
   # Illustrative shifts [L/mol], not fitted parameters.
   model = ConstantVolumeTranslation(base, [0.01, 0.02])
   volume = model.volume([0.8, 1.2], 1.0, 300.0, root="vapor")

Alternatively, calculate shifts without constructing a model using
``ConstantVolumeTranslation.psrk_correlation`` for SRK or
``ConstantVolumeTranslation.vtpr_correlation`` for PR76. Both accept
critical temperatures [K], critical pressures [bar], and experimental
critical molar volumes [L/mol], in component order:

.. code-block:: python

   tc = [190.564, 305.32]  # K
   pc = [45.99, 48.72]    # bar
   vc = [0.0986, 0.1455]  # L/mol
   shifts = ConstantVolumeTranslation.psrk_correlation(tc, pc, vc)

The caller chooses whether and how to apply the returned shifts. The PR76
correlation alone does not construct the full VTPR model.

For :math:`C=\sum_i n_i c_i`, the physical volume at fixed pressure and
temperature is the corresponding base volume minus :math:`C`. Inherited
property methods include the consistent residual Helmholtz correction.

Linear temperature-dependent translation
----------------------------------------

``LinearVolumeTranslation`` uses one intercept and temperature slope per
component: :math:`c_i(T)=a_i+b_iT`. Supply ``a`` in L/mol and ``b`` in
L/(mol K), in the base model's component order. The model copies the base EOS.

.. code-block:: python

   from yaeos import LinearVolumeTranslation, PengRobinson76

   base = PengRobinson76([190.564, 305.32], [45.99, 48.72], [0.01142, 0.0995])
   # Illustrative parameters, not fitted.
   model = LinearVolumeTranslation(base, [0.01, 0.02], [-1e-5, -2e-5])
   volume = model.volume([0.8, 1.2], 1.0, 300.0, root="vapor")

Slopes of either sign are accepted. Parameter choice matters: temperature
dependent shifts can introduce crossing isotherms. For a pure-fluid cubic EOS
with positive base thermal expansion, nonpositive slopes are sufficient to
avoid that crossover; positive slopes must be checked over the range of
interest. See Shi and Li (2016), doi:10.1016/j.fluid.2016.09.017.

.. automodule:: yaeos.models.residual_helmholtz.volume_translation
   :members:
   :show-inheritance:
