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

        import numpy as np

        from yaeos import UNIQUAC

        b = np.array(
            [
                [0.0, -526.02, -309.64],
                [318.06, 0.0, 91.532],
                [-1325.1, -302.57, 0.0],
            ]
        )

        rs = np.array([0.92, 2.1055, 3.1878])
        qs = np.array([1.4, 1.972, 2.4])

        t = 298.15

        model = UNIQUAC(qs, rs, bij=b)

        n = np.array([2.0, 2.0, 8.0])

        gammas = np.exp(model.ln_gamma(n, t)) # [8.856, 0.860, 1.425]
    """

    def __init__(
        self, qs, rs, aij=None, bij=None, cij=None, dij=None, eij=None
    ) -> None:
        self.qs = qs
        self.rs = rs

        nc = len(qs)

        if aij is not None:
            self.aij = np.array(aij, order="F")
        else:
            self.aij = np.zeros((nc, nc), order="F")

        if bij is not None:
            self.bij = np.array(bij, order="F")
        else:
            self.bij = np.zeros((nc, nc), order="F")

        if cij is not None:
            self.cij = np.array(cij, order="F")
        else:
            self.cij = np.zeros((nc, nc), order="F")

        if dij is not None:
            self.dij = np.array(dij, order="F")
        else:
            self.dij = np.zeros((nc, nc), order="F")

        if eij is not None:
            self.eij = np.array(eij, order="F")
        else:
            self.eij = np.zeros((nc, nc), order="F")

        self.id = yaeos_c.uniquac(
            qs, rs, self.aij, self.bij, self.cij, self.dij, self.eij
        )

        self.nc = len(self.qs)

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
        dij_c = ""
        eij_c = ""

        qs_c = "qs = ["
        rs_c = "rs = ["

        for i in range(self.nc):
            if i < self.nc - 1:
                qs_c += f"{self.qs[i]}_pr, "
                rs_c += f"{self.rs[i]}_pr, "
            else:
                qs_c += f"{self.qs[i]}_pr"
                rs_c += f"{self.rs[i]}_pr"

        qs_c += "]\n"
        rs_c += "]\n\n"

        for i in range(self.nc):
            aij_c += f"aij({i + 1}, :) = ["
            bij_c += f"bij({i + 1}, :) = ["
            cij_c += f"cij({i + 1}, :) = ["
            dij_c += f"dij({i + 1}, :) = ["
            eij_c += f"eij({i + 1}, :) = ["

            for j in range(self.nc):
                if j < self.nc - 1:
                    aij_c += f"{self.aij[i, j]}_pr, "
                    bij_c += f"{self.bij[i, j]}_pr, "
                    cij_c += f"{self.cij[i, j]}_pr, "
                    dij_c += f"{self.dij[i, j]}_pr, "
                    eij_c += f"{self.eij[i, j]}_pr, "
                else:
                    aij_c += f"{self.aij[i, j]}_pr"
                    bij_c += f"{self.bij[i, j]}_pr"
                    cij_c += f"{self.cij[i, j]}_pr"
                    dij_c += f"{self.dij[i, j]}_pr"
                    eij_c += f"{self.eij[i, j]}_pr"

            aij_c += "]\n"
            bij_c += "]\n"
            cij_c += "]\n"
            dij_c += "]\n"
            eij_c += "]\n"

        fcode += qs_c
        fcode += rs_c

        fcode += aij_c + "\n"
        fcode += bij_c + "\n"
        fcode += cij_c + "\n"
        fcode += dij_c + "\n"
        fcode += eij_c + "\n"

        fcode += "ge_model = setup_uniquac(qs, rs, aij, bij, cij, dij, eij)\n"
        fcode += "\n"

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
            "type(UNIQUAC) :: ge_model\n"
            "\n"
            "real(pr) :: qs(nc), rs(nc)\n"
            "real(pr) :: aij(nc,nc), bij(nc,nc), cij(nc,nc)\n"
            "real(pr) :: dij(nc,nc), eij(nc,nc)\n\n"
        )

        return fcode
