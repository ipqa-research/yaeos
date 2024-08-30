"""Cubic EoS implementations module."""

from yaeos.core import ArModel
from yaeos.lib import yaeos_c
from yaeos.models.residual_helmholtz.cubic_eos.mixing_rules import CubicMixRule


class CubicEoS(ArModel):
    """Cubic equation of state base class.

    Parameters
    ----------
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector

    Attributes
    ----------
    nc : int
        Number of components
    tc : array_like
        Critical temperatures vector [K]
    pc : array_like
        Critical pressures vector [bar]
    w : array_like
        Acentric factors vector
    """

    def __init__(
        self,
        critical_temperatures,
        critical_pressures,
        acentric_factors,
    ) -> None:
        nc = len(critical_temperatures)
        self.nc = nc
        self.tc = critical_temperatures
        self.pc = critical_pressures
        self.w = acentric_factors

    def set_mixrule(self, mixrule: CubicMixRule) -> None:
        """Set the mixing rule for the EoS.

        Parameters
        ----------
        mixrule : CubicMixRule
            Mixing rule object
        """
        self.mixrule = mixrule
        self.mixrule.set_mixrule(self.id)


class PengRobinson76(CubicEoS):
    """Peng-Robinson 1976 cubic equation of state.

    Parameters
    ----------
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    mixrule : CubicMixRule, optional
        Mixing rule object. If no provided the quadratric mixing rule (QMR)
        with zero for kij and lij parameters is set, by default None

    Attributes
    ----------
    nc : int
        Number of components
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    id : int
        EoS identifier
    mixrule : CubicMixRule
        Mixing rule object

    Example
    -------
    .. code-block:: python

        from yaeos import PengRobinson76

        tc = [190.56, 305.32]   # Critical temperatures [K]
        pc = [45.99, 48.72]     # Critical pressures [bar]
        w = [0.0115, 0.0985]    # Acentric factors

        pr76 = PengRobinson76(tc, pc, w)
    """

    name = "PengRobinson76"

    def __init__(
        self,
        critical_temperatures,
        critical_pressures,
        acentric_factors,
        mixrule: CubicMixRule = None,
    ) -> None:

        super(PengRobinson76, self).__init__(
            critical_temperatures,
            critical_pressures,
            acentric_factors,
        )
        self.id = yaeos_c.pr76(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)


class PengRobinson78(CubicEoS):
    """Peng-Robinson 1978 cubic equation of state.

    Parameters
    ----------
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    mixrule : CubicMixRule, optional
        Mixing rule object. If no provided the quadratric mixing rule (QMR)
        with zero for kij and lij parameters is set, by default None

    Attributes
    ----------
    nc : int
        Number of components
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    id : int
        EoS identifier
    mixrule : CubicMixRule
        Mixing rule object

    Example
    -------
    .. code-block:: python

        from yaeos import PengRobinson78

        tc = [190.56, 305.32]   # Critical temperatures [K]
        pc = [45.99, 48.72]     # Critical pressures [bar]
        w = [0.0115, 0.0985]    # Acentric factors

        pr78 = PengRobinson78(tc, pc, w)
    """

    name = "PengRobinson78"

    def __init__(
        self,
        critical_temperatures,
        critical_pressures,
        acentric_factors,
        mixrule: CubicMixRule = None,
    ) -> None:

        super(PengRobinson78, self).__init__(
            critical_temperatures,
            critical_pressures,
            acentric_factors,
        )
        self.id = yaeos_c.pr78(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)


class SoaveRedlichKwong(CubicEoS):
    """Soave-Redlich-Kwong cubic equation of state.

    Parameters
    ----------
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    mixrule : CubicMixRule, optional
        Mixing rule object. If no provided the quadratric mixing rule (QMR)
        with zero for kij and lij parameters is set, by default None

    Attributes
    ----------
    nc : int
        Number of components
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    id : int
        EoS identifier
    mixrule : CubicMixRule
        Mixing rule object

    Example
    -------
    .. code-block:: python

        from yaeos import SoaveRedlichKwong

        tc = [190.56, 305.32]   # Critical temperatures [K]
        pc = [45.99, 48.72]     # Critical pressures [bar]
        w = [0.0115, 0.0985]    # Acentric factors

        srk = SoaveRedlichKwong(tc, pc, w)
    """

    name = "SoaveReldichKwong"

    def __init__(
        self,
        critical_temperatures,
        critical_pressures,
        acentric_factors,
        mixrule: CubicMixRule = None,
    ) -> None:

        super(SoaveRedlichKwong, self).__init__(
            critical_temperatures,
            critical_pressures,
            acentric_factors,
        )
        self.id = yaeos_c.srk(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)


class RKPR(CubicEoS):
    """RKPR cubic equation of state.

    Parameters
    ----------
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    critical_z : array_like
        Critical compressibility factor vector
    k : array_like, optional
        k parameter, by default None
    delta_1 : array_like, optional
        delta_1 parameter, by default None
    mixrule : CubicMixRule, optional
        Mixing rule object. If no provided the quadratric mixing rule (QMR)
        with zero for kij and lij parameters is set, by default None

    Attributes
    ----------
    nc : int
        Number of components
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    zc : array_like
        Critical compressibility factor vector
    id : int
        EoS identifier
    mixrule : CubicMixRule
        Mixing rule object

    Example
    -------
    .. code-block:: python

        from yaeos import RKPR

        tc = [190.56, 305.32]   # Critical temperatures [K]
        pc = [45.99, 48.72]     # Critical pressures [bar]
        w = [0.0115, 0.0985]    # Acentric factors
        zc = [0.27, 0.28]       # Critical compressibility factor

        rkpr = RKPR(tc, pc, w, zc)
    """

    name = "RKPR"

    def __init__(
        self,
        critical_temperatures,
        critical_pressures,
        acentric_factors,
        critical_z,
        k=None,
        delta_1=None,
        mixrule: CubicMixRule = None,
    ) -> None:

        super(RKPR, self).__init__(
            critical_temperatures,
            critical_pressures,
            acentric_factors,
        )
        self.zc = critical_z

        match (k is None, delta_1 is None):
            case (True, True):
                self.id = yaeos_c.rkpr(self.tc, self.pc, self.w, self.zc)
            case (False, True):
                self.id = yaeos_c.rkpr(self.tc, self.pc, self.w, self.zc, k=k)
            case (True, False):
                self.id = yaeos_c.rkpr(
                    self.tc, self.pc, self.w, self.zc, delta_1=delta_1
                )
            case (False, False):
                self.id = yaeos_c.rkpr(
                    self.tc, self.pc, self.w, self.zc, k=k, delta_1=delta_1
                )
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)