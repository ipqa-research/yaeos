"""Solvers.

Module that contains different solvers for specific porpouses.
"""

import numpy as np

from scipy.signal import argrelmax, argrelmin

from yaeos.core import ArModel


def binary_isofugacity_x1y1pt(x, p, t, model: ArModel):
    """Isofugacity evaluation at a given P and T."""
    y1, x1 = x

    x = np.array([x1, 1 - x1])
    y = np.array([y1, 1 - y1])

    lnphi_x = model.lnphi_pt(x, pressure=p, temperature=t, root="stable")
    lnphi_y = model.lnphi_pt(y, pressure=p, temperature=t, root="stable")

    isofug = np.log(x) + lnphi_x - (np.log(y) + lnphi_y)

    return isofug**2


def solve_pt(
    model: ArModel, pressure: float, temperature: float, kind: str
) -> tuple[float, float]:
    """Solve a point at a given P and T.

    This function solves a point at a given pressure and temperature. That is,
    solve the composition of two phases in equilibrium at T and P. The function
    handles liquid-vapor or liquid-liquid.

    Parameters
    ----------
    model : ArModel
        yaeos ArModel object.
    pressure : float
        Pressure [bar].
    temperature : float
        Temperature [K].
    kind : str
        Kind of phase equilibrium calculation. Options are:
        - "liquid-liquid": liquid-liquid equilibrium calculation.
        - "PT": pressure-temperature calculation.

    Returns
    -------
    tuple[float, float]
        x1, y1: Mole fractions of component 1 (light component) in both phases.
    """
    p, t = pressure, temperature

    try:
        x10, y10 = find_init_binary_ll(model, p, t, kind)
    except ValueError:
        x10, y10 = 0.1, 0.9

    mean = (x10 + y10) / 2

    z = [mean, 1 - mean]
    y0 = np.array([y10, 1 - y10])
    x0 = np.array([x10, 1 - x10])

    flash = model.flash_pt(z, pressure=p, temperature=t, k0=y0 / x0)

    x1 = flash["x"][0]
    y1 = flash["y"][0]

    return x1, y1


def find_init_binary_ll(
    model: ArModel, pressure: float, temperature: float, kind: str
) -> tuple[float, float]:
    """Find initial guess for a binary liquid-liquid system.

    Parameters
    ----------
    model : ArModel
        yaeos ArModel object.
    pressure : float
        Pressure [bar].
    temperature : float
        Temperature [K].
    kind : str
        Kind of phase equilibrium calculation. Options are:
        - "liquid-liquid": liquid-liquid equilibrium calculation.
        - "PT": pressure-temperature calculation.

    Returns
    -------
    tuple[float, float]
        x1, y1: Initial guess for mole fractions of component 1 (light
        component) in both phases.
    """
    p, t = pressure, temperature

    if kind == "liquid-liquid":
        root = "liquid"
    else:
        root = "stable"

    zs = np.linspace(1e-15, 1 - 1e-15, 100)

    phis = np.array(
        [
            model.lnphi_pt([z, 1 - z], temperature=t, pressure=p, root=root)
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
