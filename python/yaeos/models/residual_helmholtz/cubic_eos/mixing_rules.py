"""Cubic EoS mixing rules implementations module."""

from abc import ABC, abstractmethod

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
        self.kij = kij
        self.lij = lij

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set quadratic mix rule method.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model
        """
        yaeos_c.set_qmr(ar_model_id, self.kij, self.lij)


class MHV(CubicMixRule):
    def __init__(self, ge: GeModel, q: float, lij=None) -> None:
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
        self.ge = ge
        self.q = q
        self.lij = lij

    def set_mixrule(self, ar_model_id: int) -> None:
        """Set modified Huron-Vidal mix rule method.

        Parameters
        ----------
        ar_model_id : int
            ID of the cubic EoS model
        """
        yaeos_c.set_mhv(ar_model_id, self.ge.id, self.q)
