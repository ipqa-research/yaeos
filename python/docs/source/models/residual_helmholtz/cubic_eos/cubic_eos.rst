Cubic EoS
=========

Use a concrete class such as ``PengRobinson76`` or ``SoaveRedlichKwong``;
``CubicEoS`` supplies their shared interface. Parameter arrays must have
matching lengths and component order. Critical temperatures are in K,
critical pressures in bar, and acentric factors are dimensionless.
``RKPR`` additionally needs critical compressibility factors.

For example, this builds a methane/ethane model (methane first, ethane
second) and round-trips a homogeneous state: the volume is solved at a given
pressure and temperature, and the pressure is then recomputed from that
volume. The two pressures should agree to within the solver's tolerance.

.. code-block:: python

   import numpy as np
   from yaeos import PengRobinson76

   model = PengRobinson76([190.56, 305.32], [45.99, 48.72], [0.0115, 0.0985])
   n = np.array([0.8, 0.2])  # mol; methane first, ethane second
   V = model.volume(n, 1.0, 300.0)  # L at 1 bar and 300 K
   P = model.pressure(n, V, 300.0)  # bar; should come back near 1.0
   assert np.isfinite(V) and V > 0
   assert np.isclose(P, 1.0, rtol=1e-3)

The default mixing rule for the ordinary PR/SRK constructors is quadratic
with zero interaction parameters. Zero is a starting assumption, not a
fitted result. ``PSRK`` also requires molecular groups (its
Mathias-Copeman alpha-function parameters are optional); it is not created
from critical constants alone.

Continue with :doc:`mixing_rules` and :doc:`/core` for properties, root
selection, and phase-split calculations.

.. automodule:: yaeos.models.residual_helmholtz.cubic_eos.cubic_eos
    :members:
    :undoc-members:
    :show-inheritance:
