Envelopes Module
================

An envelope follows a phase boundary, rather than repeatedly solving
unrelated flashes. Use a PT envelope to examine temperature-pressure limits
at fixed overall composition, PX at fixed temperature, and TX at fixed
pressure. Start with the
:doc:`phase-equilibrium tutorial <tutorial/equations_of_state/phase_equilibrium/phase_equilibrium>`;
the calculation methods themselves belong to :doc:`core`.

Reading an envelope
-------------------

Results are ``PTEnvelope``, ``PXEnvelope``, or ``TXEnvelope`` objects. Read
``pressures`` in bar, ``temperatures`` in K, and ``df`` for tabular inspection.
The objects also accept short string keys, which return the corresponding
array, and slice, list, or array indexing, which returns a new envelope with
the selected points (integer indices are not supported, and a TX envelope
does not accept slices). A PT envelope answers to
``"T"``, ``"P"``, ``"Tc"``, ``"Pc"``, ``"z"``, ``"x"``, and ``"w"``, while the
composition envelopes answer to ``"alpha"`` (or ``"a"``), ``"x"``, and
``"w"``, plus ``"P"`` for a PX envelope and ``"T"`` for a TX envelope.

Composition arrays retain the model's component order: the documented
``(n_points, n_phases, n_components)`` shape separates traced points, main
phases, and components. Reference-phase compositions are stored separately.
Do not flatten these axes into a single composition vector.

For composition scans, ``alphas`` parameterizes the overall feed between
``global_composition_0`` and ``global_composition_i`` as
``z = (1 - alpha) * z0 + alpha * zi``. Alpha is not generally the mole
fraction of one component, and is not the phase amount ``beta``.

Inspect iteration counts, phase kinds, and critical-point arrays alongside
the curve. A trace may stop at a numerical or user-specified limit; its last
point is not automatically a physical endpoint. The points follow the
continuation path, so a branch traced from the other end is not ordered by
increasing temperature or pressure: keep the point order when plotting.

.. automodule:: yaeos.envelopes
    :members:
    :undoc-members:
    :show-inheritance:
