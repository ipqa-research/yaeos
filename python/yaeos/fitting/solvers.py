"""Solvers

Module that contains different solvers for specific porpouses.
"""

import numpy as np


def binary_isofugacity_x1y1PT(X, P, T, model):
    """Isofugacity evaluation at a given P and T"""

    y1, x1 = X

    x = np.array([x1, 1 - x1])
    y = np.array([y1, 1 - y1])

    lnphi_x = model.lnphi_pt(x, pressure=P, temperature=T, root="stable")
    lnphi_y = model.lnphi_pt(y, pressure=P, temperature=T, root="stable")

    isofug = np.log(x) + lnphi_x - (np.log(y) + lnphi_y)

    return isofug


def solve_PT(model, P, T):
    """ """
    from scipy.optimize import root

    x10, y10 = find_init_binary_ll(model, P, T)

    X0 = [x10, y10]
    sol = root(binary_isofugacity_x1y1PT, X0, (P, T, model))
    x1, y1 = sol.x

    return x1, y1


def find_init_binary_ll(model, pressure, temperature):
    from scipy.signal import argrelmin, argrelmax

    (
        P,
        T,
    ) = (
        pressure,
        temperature,
    )

    zs = np.linspace(1e-15, 1 - 1e-15, 100)

    phis = np.array(
        [
            model.lnphi_pt(
                [z, 1 - z], temperature=T, pressure=P, root="liquid"
            )
            for z in zs
        ]
    )
    phis = np.exp(phis)
    fug_1 = zs * phis[:, 0] * P

    argmin = argrelmin(zs * phis[:, 0] * P)[-1] + 1
    argmax = argrelmax(zs * phis[:, 0] * P)[0] - 1

    fug = np.mean([fug_1[argmin], fug_1[argmax]])

    if fug > fug_1[-1]:
        fug = np.mean([fug_1[argmin[0]], fug_1[-1]])

    msk = zs < zs[argmax]
    x1 = zs[msk][np.argmin(np.abs(fug - fug_1[msk]))]

    msk = zs > zs[argmin]
    y1 = zs[msk][np.argmin(np.abs(fug - fug_1[msk]))]

    return x1, y1
