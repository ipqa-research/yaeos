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

        self.m = m
        self.sigma = sigma
        self.epsilon_k = epsilon_k
        self.kij = kij

        self.id = yaeos_c.pcsaft(m, sigma, epsilon_k, kij)

    def size(self) -> int:
        """Return the number of components in the model."""
        return len(self.m)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        # Pure component parameters
        fcode += (
            f"m = [{', '.join(str(m) + '_pr' for m in self.m)}]\n"
            f"sigma = [{', '.join(str(s) + '_pr' for s in self.sigma)}]\n"
            "epsilon_k = "
            f"[{', '.join(str(ek) + '_pr' for ek in self.epsilon_k)}]\n\n"
        )

        # Binary interaction parameters
        kij_c = ""

        for i in range(len(self.kij)):
            kij_c += f"kij({i + 1}, :) = ["

            for j in range(len(self.kij)):
                if j < len(self.kij) - 1:
                    kij_c += f"{self.kij[i][j]}_pr, "
                else:
                    kij_c += f"{self.kij[i][j]}_pr]\n"

        fcode += kij_c + "\n"

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        fcode = (
            f"integer, parameter :: nc={self.size()}\n"
            "\n"
            "type(PcSaft) :: ar_model\n"
            "\n"
            "real(pr) :: m(nc), sigma(nc), epsilon_k(nc), kij(nc,nc)\n"
            "\n"
        )

        return fcode
