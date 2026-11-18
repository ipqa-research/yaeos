Core Module
===========

.. _float: https://docs.python.org/3/library/functions.html#float

Where the shared calculations live
----------------------------------

Construct a concrete :doc:`model <models/init>` first; ``ArModel`` and
``GeModel`` describe shared interfaces rather than fluids you instantiate.
An equation-of-state object inherits ``ArModel`` methods such as ``pressure``,
``volume``, ``flash_pt``, and ``saturation_pressure``. An excess-Gibbs object
inherits ``GeModel`` methods such as ``ln_gamma`` and ``flash_t``.

Before calling a method
-----------------------

* Use the :ref:`Python input conventions <python-input-conventions>`.
  ``pressure(moles, volume, temperature)`` takes total volume in L, not
  density. ``volume(moles, pressure, temperature)`` returns total volume in L.
  Molar density is ``sum(moles) / volume`` in mol/L.
* ``lnphi_pt`` and ``ln_gamma`` return natural logarithms of dimensionless
  coefficients. Apply ``numpy.exp`` to obtain the coefficients themselves.
  Watch the argument order: on an equation-of-state model, ``ln_gamma``
  takes ``(moles, pressure, temperature)`` and can select a volume root,
  while on an excess-Gibbs model it takes ``(moles, temperature)``.
* Methods named ``residual`` or ``excess`` return a departure from a
  reference state, not a total property. Divide extensive energy results by
  total moles for molar values; add the appropriate reference contribution
  when a total property is needed.
* When derivative flags are enabled, property methods return
  ``(value, derivatives)`` instead of just ``value``. For example,
  ``p, derivatives = model.pressure(n, V, T, dt=True)`` gives
  ``derivatives["dt"]`` in bar/K. Unrequested entries may be ``None``.

Choosing a calculation
-----------------------

``volume(..., root="liquid")`` or ``root="vapor"`` selects a homogeneous
volume root; ``root="stable"`` returns whichever of the two roots has the
lower Gibbs energy. Even the stable root is not a phase-split calculation: a
mixture may lower its free energy by separating into different compositions.
Use stability analysis and a flash to investigate that possibility. Check
that the returned volume is positive and finite: some implementations signal
a failed solve with a negative value such as ``-1``.
For a two-phase ``ArModel.flash_pt`` result, ``x`` and ``y`` are heavy- and
light-phase mole fractions; ``beta`` is the light-phase mole fraction of the
whole feed, not its volume fraction. The documented failure value is
``beta == -1``. For an admissible split, check that the values are finite,
that the phase compositions are normalized, and that the material balance
``z = (1 - beta) * x + beta * y`` holds. A converged split alone does not
prove global stability.

``Vx`` and ``Vy`` in the flash and saturation results are phase **molar**
volumes in L/mol, as documented by the Fortran API and used in the tutorial
notebooks. They are therefore not directly comparable with the total volume
returned by ``volume(...)``, which is the volume of the whole ``moles``
vector in L. Note also that
``flash_vt`` takes a **molar** volume in L/mol, while ``pressure`` and
``volume`` use a total volume in L.

For ``GeModel.flash_t``, ``x`` and ``y`` label two liquid phases and ``beta``
is the fraction of the second phase; ``y`` is not a vapor composition.
Continue with the
:doc:`phase-equilibrium tutorial <tutorial/equations_of_state/phase_equilibrium/phase_equilibrium>`
or the :doc:`liquid-mixture tutorials <tutorial/gibbs_excess/gibbs_excess>`.

.. automodule:: yaeos.core
    :members:
    :undoc-members:
    :show-inheritance:
