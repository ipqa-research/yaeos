"""Solvers.

Module that contains different solvers for specific porpouses.
"""

import numpy as np


def binary_isofugacity_x1y1pt(x, p, t, model):
    """Isofugacity evaluation at a given P and T."""

    y1, x1 = x

    x = np.array([x1, 1 - x1])
    y = np.array([y1, 1 - y1])

    lnphi_x = model.lnphi_pt(x, pressure=p, temperature=t, root="stable")
    lnphi_y = model.lnphi_pt(y, pressure=p, temperature=t, root="stable")

    isofug = np.log(x) + lnphi_x - (np.log(y) + lnphi_y)

    return isofug


def solve_pt(model, p, t):
    """Solve a point at a given P and T."""
    from scipy.optimize import root

    x10, y10 = find_init_binary_ll(model, p, t)

    x0 = [x10, y10]
    sol = root(binary_isofugacity_x1y1pt, x0, (p, t, model))
    x1, y1 = sol.x

    return x1, y1


def find_init_binary_ll(model, pressure, temperature):
    """Find initial guess for a binary liquid-liquid system."""
    from scipy.signal import argrelmin, argrelmax

    (
        p,
        t,
    ) = (
        pressure,
        temperature,
    )

    zs = np.linspace(1e-15, 1 - 1e-15, 100)

    phis = np.array(
        [
            model.lnphi_pt(
                [z, 1 - z], temperature=t, pressure=p, root="liquid"
            )
            for z in zs
        ]
    )
    phis = np.exp(phis)
    fug_1 = zs * phis[:, 0] * p

    argmin = argrelmin(zs * phis[:, 0] * p)[-1] + 1
    argmax = argrelmax(zs * phis[:, 0] * p)[0] - 1

    fug = np.mean([fug_1[argmin], fug_1[argmax]])

    if fug > fug_1[-1]:
        fug = np.mean([fug_1[argmin[0]], fug_1[-1]])

    msk = zs < zs[argmax]
    x1 = zs[msk][np.argmin(np.abs(fug - fug_1[msk]))]

    msk = zs > zs[argmin]
    y1 = zs[msk][np.argmin(np.abs(fug - fug_1[msk]))]

    return x1, y1
