"""Non-random two-liquid model (NRTL) module."""

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
        self.a = a
        self.b = b
        self.c = c
        self.id = yaeos_c.nrtl(a, b, c)

        self.nc = len(a)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        aij_c = ""
        bij_c = ""
        cij_c = ""

        for i in range(self.nc):
            aij_c += f"aij({i + 1}, :) = ["
            bij_c += f"bij({i + 1}, :) = ["
            cij_c += f"cij({i + 1}, :) = ["

            for j in range(self.nc):
                if j < self.nc - 1:
                    aij_c += f"{self.a[i, j]}_pr, "
                    bij_c += f"{self.b[i, j]}_pr, "
                    cij_c += f"{self.c[i, j]}_pr, "
                else:
                    aij_c += f"{self.a[i, j]}_pr"
                    bij_c += f"{self.b[i, j]}_pr"
                    cij_c += f"{self.c[i, j]}_pr"

            aij_c += "]\n"
            bij_c += "]\n"
            cij_c += "]\n"

        fcode += aij_c + "\n"
        fcode += bij_c + "\n"
        fcode += cij_c + "\n"

        fcode += "ge_model = NRTL(aij, bij, cij)\n\n"

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        fcode = (
            f"integer, parameter :: nc={self.nc}\n"
            "\n"
            "type(NRTL) :: ge_model\n"
            "\n"
            "real(pr) :: aij(nc,nc), bij(nc,nc), cij(nc,nc)\n\n"
        )

        return fcode
