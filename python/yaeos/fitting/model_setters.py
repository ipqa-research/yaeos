"""Model Setters.

Compilation of functions that set a model's parameters to the values
given in the input arguments.
"""

import numpy as np


def fit_kij_lij(x, model, fit_kij, fit_lij):
    """Set the kij and/or lij parameter of the model.

    Intended to be used as `model_setter` in
    :class:`~yaeos.fitting.core.BinaryFitter` with
    ``model_setter_args=(model, fit_kij, fit_lij)``.

    Parameters
    ----------
    x : array_like
        Optimization parameters: ``[kij, lij]`` if both are fitted,
        otherwise ``[kij]`` or ``[lij]``.
    model : CubicEoS
        Binary cubic EoS whose mixing rule has `kij` and `lij` attributes
        (e.g. :class:`~yaeos.models.residual_helmholtz.cubic_eos.QMR`).
    fit_kij : bool
        Whether to fit the kij parameter.
    fit_lij : bool
        Whether to fit the lij parameter.

    Returns
    -------
    CubicEoS
        The same model object with its mixing rule updated.
    """
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
    """Set an MHV mixing rule with NRTL GE in a binary cubic EoS.

    Parameters
    ----------
    x : array_like
        Optimization parameters ``[a12, a21, b12, b21, alpha]`` of the NRTL
        model.
    args : tuple
        ``(model, q)``: the cubic EoS and the MHV `q` parameter. Note that
        :class:`~yaeos.fitting.core.BinaryFitter` calls the setter as
        ``model_setter(x, *model_setter_args)``, so this tuple must be
        wrapped: ``model_setter_args=((model, q),)``.

    Returns
    -------
    CubicEoS
        The same model object with its mixing rule updated.
    """
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
    """Set an HV mixing rule with NRTL GE in a binary cubic EoS.

    Parameters
    ----------
    x : array_like
        Optimization parameters ``[a12, a21, b12, b21, alpha]`` of the NRTL
        model.
    args : tuple
        ``(model,)``: the cubic EoS. Note that
        :class:`~yaeos.fitting.core.BinaryFitter` calls the setter as
        ``model_setter(x, *model_setter_args)``, so this tuple must be
        wrapped: ``model_setter_args=((model,),)``.

    Returns
    -------
    CubicEoS
        The same model object with its mixing rule updated.
    """
    from yaeos.models.excess_gibbs import NRTL
    from yaeos.models.residual_helmholtz.cubic_eos.mixing_rules import HV

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
