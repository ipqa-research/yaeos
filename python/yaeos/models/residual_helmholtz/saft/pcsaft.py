"""PC-SAFT Equation of State."""

import numpy as np

from yaeos.core import ArModel
from yaeos.lib import yaeos_c


class PCSAFT(ArModel):
    """PC-SAFT Equation of State.

    This class implements the residual contribution of the PC-SAFT equation of
    state for multi-component systems.

    Parameters
    ----------
    m: list, float
        Segment number for each component.
    sigma: list, float
        Segment diameter for each component (in Angstroms).
    epsilon_k: list, float
        Segment energy parameter for each component (in Kelvin).
    kij: list, list, float, optional
        Binary interaction parameter matrix. Default is None, which sets all
        interaction parameters to zero.

    Example
    -------
    .. code-block:: python
        from yaeos import PCSAFT

        m = [1.0582, 3.3004]
        sigma = [3.6316, 3.8639]
        epsilon_k = [145.5257, 224.0780]

        model = PCSAFT(m, sigma, epsilon_k)

        # Or use kij matrix for binary interaction parameters
        kij = [[0.0, 0.03],
               [0.03, 0.0]]
        model = PCSAFT(m, sigma, epsilon_k, kij=kij)
    """

    def __init__(
        self, m: np.ndarray, sigma: np.ndarray, epsilon_k: np.ndarray, kij=None
    ):
        """Initialize PC-SAFT model."""
        if kij is None:
            kij = [[0.0 for _ in m] for _ in m]
        self.size = len(m)
        self.id = yaeos_c.pcsaft(m, sigma, epsilon_k, kij)

    def size(self) -> int:
        """Return the number of components in the model."""
        return self.size
