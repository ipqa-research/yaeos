"""UNIQUAC (UNIversal QUAsiChemical) Excess Gibbs free energy model."""

import numpy as np

from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class UNIQUAC(GeModel):
    """UNIQUAC (UNIversal QUAsiChemical) Excess Gibbs free energy model.

    Please refer to the `yaeos` user documentation for an in-depth look at the
    model's information: https://ipqa-research.github.io/yaeos/page/index.html

    Parameters
    ----------
    qs : array_like
        Molecule's relative areas :math:`Q_i`
    rs : array_like
        Molecule's relative volumes :math:`R_i`
    aij : array_like
        Interaction parameters matrix :math:`a_{ij}` zero matrix if no
        provided, by default None
    bij : array_like
        Interaction parameters matrix :math:`b_{ij}` zero matrix if no
        provided, by default None
    cij : array_like
        Interaction parameters matrix :math:`c_{ij}` zero matrix if no
        provided, by default None
    dij : array_like
        Interaction parameters matrix :math:`d_{ij}` zero matrix if no
        provided, by default None
    eij : array_like
        Interaction parameters matrix :math:`e_{ij}` zero matrix if no
        provided, by default None

    Attributes
    ----------
    qs : array_like
        Molecule's relative areas :math:`Q_i`
    rs : array_like
        Molecule's relative volumes :math:`R_i`
    aij : array_like
        Interaction parameters matrix :math:`a_{ij}`
    bij : array_like
        Interaction parameters matrix :math:`b_{ij}`
    cij : array_like
        Interaction parameters matrix :math:`c_{ij}`
    dij : array_like
        Interaction parameters matrix :math:`d_{ij}`
    eij : array_like
        Interaction parameters matrix :math:`e_{ij}`

    Example
    -------
    .. code-block:: python

        from yaeos import UNIFACVLE

        # Groups for water and ethanol
        water = {16: 1}
        ethanol = {1: 1, 2: 1, 14: 1}

        groups = [water, ethanol]

        model = UNIFAVLE(groups)

        model.ln_gamma([0.5, 0.5], 298.15)
    """

    def __init__(
        self, qs, rs, aij=None, bij=None, cij=None, dij=None, eij=None
    ) -> None:
        self.qs = qs
        self.rs = rs

        nc = len(qs)

        if aij is not None:
            self.aij = aij
        else:
            self.aij = np.zeros((nc, nc), order="F")

        if bij is not None:
            self.bij = bij
        else:
            self.bij = np.zeros((nc, nc), order="F")

        if cij is not None:
            self.cij = cij
        else:
            self.cij = np.zeros((nc, nc), order="F")

        if dij is not None:
            self.dij = dij
        else:
            self.dij = np.zeros((nc, nc), order="F")

        if eij is not None:
            self.eij = eij
        else:
            self.eij = np.zeros((nc, nc), order="F")

        self.id = yaeos_c.uniquac(
            qs, rs, self.aij, self.bij, self.cij, self.dij, self.eij
        )
