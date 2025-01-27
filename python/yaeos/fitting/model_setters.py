"""Model Setters.

Compilation of functions that set a model's parameters to the values
given in the input arguments.
"""

import numpy as np


def fit_kij_lij(x, args):
    """Set the kij and/or lij parameter of the model."""
    model = args[0]
    fit_kij, fit_lij = args[1:]

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
