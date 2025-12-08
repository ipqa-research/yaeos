"""Non-random two-liquid model (NRTL) module."""

import numpy as np

from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class NRTL(GeModel):
    """Non-random two-liquid model (NRTL) class.

    Please refer to the `yaeos` user documentation for an in-depth look at the
    model's information: https://ipqa-research.github.io/yaeos/page/index.html

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
        self.a = np.array(a, order="F")
        self.b = np.array(b, order="F")
        self.c = np.array(c, order="F")
        self.id = yaeos_c.nrtl(a, b, c)

    def size(self) -> int:
        """Get the number of components.

        Returns
        -------
        int
            Number of components
        """
        return self.a.shape[0]
