Constant volume translation
===========================

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

.. automodule:: yaeos.models.residual_helmholtz.volume_translation
   :members:
   :show-inheritance:
