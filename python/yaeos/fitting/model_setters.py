"""Model Setters.

Compilation of functions that set a model's parameters to the values
given in the input arguments.
"""

import numpy as np


def fit_kij_lij(x, model, fit_kij, fit_lij):
    """Set the kij and/or lij parameter of the model."""
    mr = model.mixrule

    if fit_kij and fit_lij:
        kij = x[0]
        lij = x[1]
        mr.kij = np.array([[0, kij], [kij, 0]], order="F")
        mr.lij = np.array([[0, lij], [lij, 0]], order="F")
    elif fit_kij:
        kij = x[0]
        mr.kij = np.array([[0, kij], [kij, 0]], order="F")
    elif fit_lij:
        lij = x[0]
        mr.lij = np.array([[0, lij], [lij, 0]], order="F")

    model.set_mixrule(mr)
    return model


def fit_mhv_nrtl(x, args):
    """Fit the MHV mixing rule for Cubic EoS with NRTL GE."""
    from yaeos.models import NRTL, MHV

    a12, a21, b12, b21, alpha = x

    model = args[0]
    q = args[1]

    a = [
        [
            0,
            a12,
        ],
        [a21, 0],
    ]
    b = [[0, b12], [b21, 0]]
    c = [[0, alpha], [alpha, 0]]

    nrtl = NRTL(a, b, c)
    mr = MHV(ge=nrtl, q=q)
    model.set_mixrule(mr)
    return model


def fit_hv_nrtl(x, args):
    """Fit the HV mixing rule for Cubic EoS with NRTL GE."""
    from yaeos.models import NRTL, HV

    a12, a21, b12, b21, alpha = x

    model = args[0]

    a = [
        [
            0,
            a12,
        ],
        [a21, 0],
    ]
    b = [[0, b12], [b21, 0]]
    c = [[0, alpha], [alpha, 0]]

    nrtl = NRTL(a, b, c)
    mr = HV(nrtl)
    model.set_mixrule(mr)
    return model
