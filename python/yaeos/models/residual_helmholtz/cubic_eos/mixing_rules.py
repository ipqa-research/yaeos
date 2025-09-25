"""Cubic EoS mixing rules implementations module."""

from abc import ABC, abstractmethod

import numpy as np

from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class CubicMixRule(ABC):
    """Cubic mix rule abstract class."""

    @abstractmethod
    def set_mixrule(self, ar_model_id: int) -> None:
        """Set mix rule abstract method.

        Changes the default mixing rule of the given cubic EoS model.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model

        Raises
        ------
        NotImplementedError
            Abstract error, this method must be implemented in the subclass
        """
        raise NotImplementedError

    @abstractmethod
    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        pass

    @abstractmethod
    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        pass


class QMR(CubicMixRule):
    """Quadratic mixing rule.

    Parameters
    ----------
    kij : array_like
        kij binary interaction parameters matrix
    lij : array_like
        lij binary interaction parameters matrix

    Attributes
    ----------
    kij : array_like
        kij binary interaction parameters matrix
    lij : array_like
        lij binary interaction parameters matrix

    Example
    -------
    .. code-block:: python

        from yaeos import QMR, SoaveRedlichKwong

        kij = [[0.0, 0.1], [0.1, 0.0]]
        lij = [[0.0, 0.02], [0.02, 0.0]]

        mixrule = QMR(kij, lij)     # Quadratic mixing rule instance

        tc = [305.32, 469.7]        # critical temperature [K]
        pc = [48.72, 33.7]          # critical pressure [bar]
        w = [0.0995, 0.152]         # acentric factor

        model = SoaveRedlichKwong(tc, pc, w, mixrule)
    """

    def __init__(self, kij, lij) -> None:
        self.kij = np.array(kij, order="F")
        self.lij = np.array(lij, order="F")

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set quadratic mix rule method.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model
        """
        yaeos_c.set_qmr(ar_model_id, self.kij, self.lij)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        kij_c = ""
        lij_c = ""

        for i in range(len(self.kij)):
            kij_c += f"kij({i + 1}, :) = ["
            lij_c += f"lij({i + 1}, :) = ["

            for j in range(len(self.kij)):
                if j < len(self.kij) - 1:
                    kij_c += f"{self.kij[i][j]}_pr, "
                    lij_c += f"{self.lij[i][j]}_pr, "
                else:
                    kij_c += f"{self.kij[i][j]}_pr]\n"
                    lij_c += f"{self.lij[i][j]}_pr]\n"

        fcode += kij_c + "\n"
        fcode += lij_c + "\n"
        fcode += "mixrule = QMR(k=kij, l=lij)\n\n"

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """

        fcode = (
            "type(QMR) :: mixrule"
            "\n"
            "real(pr) :: kij(nc, nc), lij(nc, nc)\n\n"
        )

        return fcode


class QMRTD(CubicMixRule):
    r"""Quadratic mixing rule, with temperature dependence.

    ..math::
        k_{ij}(T) = k_{ij}^{\infty} + k_{ij}^0 \exp{\left(-T/T^{ref}\right)}

    Parameters
    ----------
    kij_0 : array_like
        kij_0 binary interaction parameters matrix
    kij_inf : array_like
        kij_inf binary interaction parameters matrix
    t_ref: array_like
        Reference temperature
    lij : array_like
        lij binary interaction parameters matrix

    Attributes
    ----------
    kij_0 : array_like
        kij_0 binary interaction parameters matrix
    kij_inf : array_like
        kij_inf binary interaction parameters matrix
    t_ref: array_like
        Reference temperature
    lij : array_like
        lij binary interaction parameters matrix

    Example
    -------
    .. code-block:: python

        from yaeos import QMRTD, SoaveRedlichKwong

        kij_0 = [[0.0, 0.1], [0.1, 0.0]]
        kij_inf = [[0.0, 0.1], [0.1, 0.0]]
        Tref = [[0.0, 390], [390, 0.0]]
        lij = [[0.0, 0.02], [0.02, 0.0]]

        # Quadratic mixing rule instance
        mixrule = QMRTD(kij_0, kij_inf, Tref, lij)

        tc = [305.32, 469.7]        # critical temperature [K]
        pc = [48.72, 33.7]          # critical pressure [bar]
        w = [0.0995, 0.152]         # acentric factor

        model = SoaveRedlichKwong(tc, pc, w, mixrule)
    """

    def __init__(self, kij_0, kij_inf, t_ref, lij) -> None:
        self.kij_0 = np.array(kij_0, order="F")
        self.kij_inf = np.array(kij_inf, order="F")
        self.t_ref = np.array(t_ref, order="F")
        self.lij = np.array(lij, order="F")

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set mix rule method."""
        yaeos_c.set_qmrtd(
            ar_model_id,
            kij_0=self.kij_0,
            kij_inf=self.kij_inf,
            t_star=self.t_ref,
            lij=self.lij,
        )

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        kij0_c = ""
        kijinf_c = ""
        lij_c = ""
        t_ref_c = ""

        for i in range(len(self.kij_0)):
            kij0_c += f"kij_0({i + 1}, :) = ["
            kijinf_c += f"kij_inf({i + 1}, :) = ["
            lij_c += f"lij({i + 1}, :) = ["
            t_ref_c += f"t_ref({i + 1}, :) = ["

            for j in range(len(self.kij_0)):
                if j < len(self.kij_0) - 1:
                    kij0_c += f"{self.kij_0[i][j]}_pr, "
                    kijinf_c += f"{self.kij_inf[i][j]}_pr, "
                    lij_c += f"{self.lij[i][j]}_pr, "
                    t_ref_c += f"{self.t_ref[i][j]}_pr, "
                else:
                    kij0_c += f"{self.kij_0[i][j]}_pr]\n"
                    kijinf_c += f"{self.kij_inf[i][j]}_pr]\n"
                    lij_c += f"{self.lij[i][j]}_pr]\n"
                    t_ref_c += f"{self.t_ref[i][j]}_pr]\n"

        fcode += kij0_c + "\n"
        fcode += kijinf_c + "\n"
        fcode += lij_c + "\n"
        fcode += t_ref_c + "\n"

        fcode += "mixrule = QMRTD(k=kij_inf, k0=kij_0, Tref=t_ref, l=lij)\n\n"

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """

        fcode = (
            "type(QMRTD) :: mixrule"
            "\n"
            "real(pr) :: kij_inf(nc, nc), kij_0(nc, nc), lij(nc, nc)\n"
            "real(pr) :: t_ref(nc, nc)\n\n"
        )

        return fcode


class MHV(CubicMixRule):
    """Modified Huron-Vidal mixing rule.

    Parameters
    ----------
    ge : GeModel
        Excess Gibbs energy model
    q : float
        q parameter. Use:
            q = -0.594 for Soave-Redlich-Kwong
            q = -0.53 for Peng-Robinson
            q = -0.85 for Van der Waals
    lij : array_like, optional
        lij binary interaction parameters matrix, by default None

    Attributes
    ----------
    ge : GeModel
        Excess Gibbs energy model
    q : float
        q parameter
    lij : array_like
        lij binary interaction parameters matrix

    Example
    -------
    .. code-block:: python

        from yaeos import MHV, SoaveRedlichKwong, NRTL

        tc = [647.14, 513.92]               # critical temperature [K]
        pc = [220.64, 61.48]                # critical pressure [bar]
        w =  [0.344, 0.649]                 # acentric factor

        a = [[0, 3.458], [-0.801, 0]]       # NRTL aij parameters
        b = [[0, -586.1], [246.2, 0]]       # NRTL bij parameters
        c = [[0, 0.3], [0.3, 0]]            # NRTL cij parameters

        ge_model = NRTL(a, b, c)
        mixrule = MHV(ge_model, q=-0.53)

        model_mhv = PengRobinson76(tc, pc, w, mixrule)
    """

    def __init__(self, ge: GeModel, q: float, lij=None) -> None:
        self.ge = ge
        self.q = q
        self.lij = np.array(lij, order="F")

        self._have_lij = True if lij is not None else False

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set modified Huron-Vidal mix rule method.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model
        """
        yaeos_c.set_mhv(ar_model_id, self.ge.id, self.q)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        fcode += self.ge._model_params_as_str()

        lij_c = ""

        if self._have_lij:
            for i in range(self.ge.nc):
                lij_c += f"lij({i + 1}, :) = ["

                for j in range(self.ge.nc):
                    if j < self.ge.nc - 1:
                        lij_c += f"{self.lij[i][j]}_pr, "
                    else:
                        lij_c += f"{self.lij[i][j]}_pr]\n"

            fcode += lij_c + "\n"

        if self._have_lij:
            # TODO: include lij in MHV
            fcode += (
                f"mixrule = MHV(ge=ge_model, q={self.q}_pr, b=ar_model%b)\n\n"
            )
        else:
            fcode += (
                f"mixrule = MHV(ge=ge_model, q={self.q}_pr, bi=ar_model%b)\n\n"
            )

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        fcode = ""

        fcode += self.ge._model_params_declaration_as_str()

        # Just in case all possible replaces
        fcode = fcode.replace(f"integer, parameter :: nc={self.ge.nc}\n\n", "")
        fcode = fcode.replace(f"integer, parameter :: nc={self.ge.nc}\n", "")
        fcode = fcode.replace(f"integer, parameter :: nc={self.ge.nc}", "")

        fcode += "type(MHV) :: mixrule" "\n"

        if self._have_lij:
            fcode += "real(pr) :: lij(nc, nc)\n\n"

        return fcode


class HV(CubicMixRule):
    """Huron-Vidal mixing rule.

    Parameters
    ----------
    ge : GeModel
        Excess Gibbs energy model

    Attributes
    ----------
    ge : GeModel
        Excess Gibbs energy model

    Example
    -------
    .. code-block:: python

        from yaeos import HV, SoaveRedlichKwong, NRTL

        tc = [647.14, 513.92]               # critical temperature [K]
        pc = [220.64, 61.48]                # critical pressure [bar]
        w =  [0.344, 0.649]                 # acentric factor

        a = [[0, 3.458], [-0.801, 0]]       # NRTL aij parameters
        b = [[0, -586.1], [246.2, 0]]       # NRTL bij parameters
        c = [[0, 0.3], [0.3, 0]]            # NRTL cij parameters

        ge_model = NRTL(a, b, c)
        mixrule = HV(ge_model)

        model_hv = PengRobinson76(tc, pc, w, mixrule)
    """

    def __init__(self, ge: GeModel) -> None:
        self.ge = ge

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set modified Huron-Vidal mix rule method.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model
        """
        yaeos_c.set_hv(ar_model_id, self.ge.id)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        fcode += self.ge._model_params_as_str()

        fcode += (
            "mixrule = HV(ge=ge_model, "
            "bi=ar_model%b, del1=ar_model%del1)\n\n"
        )

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        fcode = ""

        fcode += self.ge._model_params_declaration_as_str()

        # Just in case all possible replaces
        fcode = fcode.replace(f"integer, parameter :: nc={self.ge.nc}\n\n", "")
        fcode = fcode.replace(f"integer, parameter :: nc={self.ge.nc}\n", "")
        fcode = fcode.replace(f"integer, parameter :: nc={self.ge.nc}", "")

        fcode += "type(HV) :: mixrule" "\n"

        return fcode


class HVNRTL(CubicMixRule):
    """Huron-Vidal mixing rule using an NRTL model.

    Huron Vidal mixing rule coupled with the excess Gibbs energy model defined
    by Huron-Vidal. This model can use the classic VdW1f mixing rules kij
    parameters when desired.

    Parameters
    ----------
    alpha : matrix_like
        NRTL alpha parameters matrix
    gji : matrix_like
        NRTL gij parameters matrix
    use_kij : matrix_like
        Boolean matrix indicating whether to use kij parameters
    kij : matrix_like
        kij binary interaction parameters matrix

    Attributes
    ----------
    alpha : matrix_like
        NRTL alpha parameters matrix
    gji : matrix_like
        NRTL gij parameters matrix
    use_kij : matrix_like
        Boolean matrix indicating whether to use kij parameters
    kij : matrix_like
        kij binary interaction parameters matrix
    """

    def __init__(self, alpha, gji, use_kij, kij) -> None:
        self.alpha = alpha
        self.gji = gji
        self.use_kij = np.array(use_kij, order="F")
        self.kij = np.array(kij, order="F")

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set modified Huron-Vidal mix rule method.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model
        """
        yaeos_c.set_hvnrtl(
            ar_model_id, self.alpha, self.gji, self.use_kij, self.kij
        )

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        # TODO
        return ""

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        # TODO
        return ""
