"""Non-random two-liquid model (NRTL) module."""

from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class NRTL(GeModel):
    """Non-random two-liquid model (NRTL) class.

    Parameters
    ----------
    a : array_like
        NRTL aij parameters matrix
    b : array_like
        NRTL bij parameters matrix
    c : array_like
        NRTL cij parameters matrix

    Attributes
    ----------
    a : array_like
        NRTL aij parameters matrix
    b : array_like
        NRTL bij parameters matrix
    c : array_like
        NRTL cij parameters matrix
    id : int
        NRTL model ID

    Example
    -------
    .. code-block:: python

        import numpy as np

        from yaeos import NRTL

        a = np.array([[0, 0.3], [0.3, 0]])
        b = np.array([[0, 0.4], [0.4, 0]])
        c = np.array([[0, 0.5], [0.5, 0]])

        nrtl = NRTL(a, b, c)
    """

    def __init__(self, a, b, c) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.id = yaeos_c.nrtl(a, b, c)
