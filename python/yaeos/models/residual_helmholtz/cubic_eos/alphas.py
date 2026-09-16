"""Alpha Functions for CubicEoS

This module contain the setter functions for the different alpha
functions existing in yaeos.
"""

from abc import ABC, abstractmethod

from yaeos.lib import yaeos_c
from yaeos.tools.writers import f1darray_as_str


class CubicAlphaFunction(ABC):
    """Cubic Alpha Function abstract class."""

    @abstractmethod
    def set_alpha(self, ar_model_id: int) -> None:
        """Set the Alpha Function abstract method.

        Changes the default Alpha Function of the given cubic EoS model.

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


class AlphaSoave(CubicAlphaFunction):
    """Soave's Alpha Function.

    Parameters
    ----------
    k: array_like
        k constant for the function

    Attributes
    ----------
    k: array_like
        k constant for the function

    Example
    -------
    .. code-block: python

        from yaeos import AlphaSoave, SoaveRedlichKwong

        k= [0.5, 0.6]

        alphasoave = AlphaSoave(k)  # AlphaSoave with custom parameters

        tc = [305.32, 469.7]        # critical temperature [K]
        pc = [48.72, 33.7]          # critical pressure [bar]
        w = [0.0995, 0.152]         # acentric factor

        model = SoaveRedlichKwong(tc, pc, w)
        model.set_alpha(alphasoave)
    """

    name = "AlphaSoave"

    def __init__(self, k):
        self.k = k

    def set_alpha(self, ar_model_id: int):
        yaeos_c.set_alpha_soave(ar_model_id, self.k)

    def _model_params_as_str(self):
        fcode = f1darray_as_str(self.k, "k")

        return fcode

    def _model_params_declaration_as_str(self):
        fcode = "type(AlphaSoave) :: alpha_function \n"
        fcode += "real(pr) :: k(nc)\n"
        return fcode


class AlphaRKPR(CubicAlphaFunction):
    """RKPR's Alpha Function.

    Parameters
    ----------
    k: array_like
        k constant for the function

    Attributes
    ----------
    k: array_like
        k constant for the function

    Example
    -------
    .. code-block: python

        from yaeos import AlphaRKPR, SoaveRedlichKwong

        k= [0.5, 0.6]

        alphaRKPR = AlphaRKPR(k)    # AlphaRKPR with custom parameters

        tc = [305.32, 469.7]        # critical temperature [K]
        pc = [48.72, 33.7]          # critical pressure [bar]
        w = [0.0995, 0.152]         # acentric factor

        model = SoaveRedlichKwong(tc, pc, w)
        model.set_alpha(alphaRKPR)
    """

    name = "AlphaRKPR"

    def __init__(self, k):
        self.k = k

    def set_alpha(self, ar_model_id: int):
        yaeos_c.set_alpha_rkpr(ar_model_id, self.k)

    def _model_params_as_str(self):
        fcode = f1darray_as_str(self.k, "k")

        return fcode

    def _model_params_declaration_as_str(self):
        fcode = "type(AlphaRKPR) :: alpha_function \n"
        fcode += "real(pr) :: k(nc)\n"
        return fcode


class AlphaMathiasCopeman(CubicAlphaFunction):
    """Mathias Copeman's Alpha Function.

    Parameters
    ----------
    c1: array_like
        c1 constant for the function
    c2: array_like
        c2 constant for the function
    c3: array_like
        c3 constant for the function

    Attributes
    ----------
    c1: array_like
        c1 constant for the function
    c2: array_like
        c2 constant for the function
    c3: array_like
        c3 constant for the function

    Example
    -------
    .. code-block: python

        from yaeos import AlphaRKPR, SoaveRedlichKwong

        c1 = [0.5, 0.6]
        c2 = [3.5, 1.6]
        c3 = [0.0, 0.0]

        # AlphaMathiasCopeman with custom parameters
        alphamc = AlphaMathiasCopeman(k)

        tc = [305.32, 469.7]        # critical temperature [K]
        pc = [48.72, 33.7]          # critical pressure [bar]
        w = [0.0995, 0.152]         # acentric factor

        model = SoaveRedlichKwong(tc, pc, w)
        model.set_alpha(alphamc)
    """

    name = "AlphaMathiasCopeman"

    def __init__(self, c1, c2, c3):
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3

    def set_alpha(self, ar_model_id: int):
        yaeos_c.set_alpha_mathiascopeman(
            ar_model_id, self.c1, self.c2, self.c3
        )

    def _model_params_as_str(self):
        c1 = f1darray_as_str(self.c1, "c1")
        c2 = f1darray_as_str(self.c2, "c2")
        c3 = f1darray_as_str(self.c3, "c3")

        fcode = "\n".join([c1, c2, c3])

        return fcode

    def _model_params_declaration_as_str(self):
        fcode = "type(AlphaMathiasCopeman) :: alpha_function \n"
        fcode += "real(pr) :: c1(nc), c2(nc), c3(nc) \n"
        return fcode
