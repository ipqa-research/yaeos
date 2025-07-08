"""Cubic EoS implementations module."""

import numpy as np

from yaeos.core import ArModel
from yaeos.lib import yaeos_c
from yaeos.models.groups import groups_from_dicts
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
    mixrule : CubicMixRule
        Cubic Mixing rule object
    """

    name = None

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
        self.mixrule = None

    def set_mixrule(self, mixrule: CubicMixRule) -> None:
        """Set the mixing rule for the EoS.

        Parameters
        ----------
        mixrule : CubicMixRule
            Mixing rule object
        """
        self.mixrule = mixrule
        self.mixrule.set_mixrule(self.id)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = (
            f"tc = [{', '.join(str(t) + '_pr' for t in self.tc)}]\n"
            f"pc = [{', '.join(str(p) + '_pr' for p in self.pc)}]\n"
            f"w = [{', '.join(str(w) + '_pr' for w in self.w)}]\n"
            "\n"
            f"model = {self.name}(tc, pc, w)"
        )

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
            "class(ArModel), allocatable :: model\n"
            "\n"
            f"real(pr) :: tc(nc), pc(nc), w(nc)\n"
        )

        return fcode


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

        self._fname = "PengRobinson76"


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

        self._fname = "PengRobinson78"


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


class PSRK(CubicEoS):
    """Predictive-Soave-Redlich-Kwong cubic equation of state.

    Parameters
    ----------
    critical_temperatures : array_like
        Critical temperatures vector [K]
    critical_pressures : array_like
        Critical pressures vector [bar]
    acentric_factors : array_like
        Acentric factors vector
    molecules: list of dict
        List of dicts with the groups and their amounts for each molecule
    c1: array_like
        Mathias-Copeman parameters c1
    c2: array_like
        Mathias-Copeman parameters c3
    c3: array_like
        Mathias-Copeman parameters c3

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

    Example
    -------
    .. code-block:: python

        from yaeos import PSRK

        # methanol/n-hexane mixture
        tc = [512.5, 507.6]   # Critical temperatures [K]
        pc = [80.84, 30.25]     # Critical pressures [bar]
        w = [0.565831, 0.301261]    # Acentric factors

        # Methanol: 1 CH3OH subgroup
        # n-hexane: 2 CH3, 4 CH2 subgroups
        molecules = [{15:1}, {1:2, 2:4}]

        # Mathias-copeman constants
        c1 = [1.31458917, 0.93830213]
        c2 = [0.0, 0.0]
        c3 = [0.0, 0.0]

        psrk = PSRK(tc, pc, w, molecules=molecules, c1=c1, c2=c2, c3=c3)

    .. code-block:: python

        # The dictionary of molecules can be created using the `ugropy` library

        import ugropy as ug

        molecules = ["methanol", "ethanol"]

        groups = [ug.get_groups(ug.psrk, molecule) for molecule in molecules]
        molecules = [
            ug.writers.to_thermo(grp.subgroups, ug.psrk) for grp in groups
        ]

        print(molecules)

        >>> [{15: 1}, {1: 1, 2: 1, 14: 1}]
    """

    name = "PSRK"

    def __init__(
        self,
        critical_temperatures,
        critical_pressures,
        acentric_factors,
        molecules,
        c1=None,
        c2=None,
        c3=None,
    ) -> None:

        super(PSRK, self).__init__(
            critical_temperatures,
            critical_pressures,
            acentric_factors,
        )
        (number_of_groups, groups_ids, groups_ammounts) = groups_from_dicts(
            molecules
        )

        if c1 is None:
            self.w = np.array(self.w)
            c1 = 0.48 + 1.574 * self.w - 0.175 * self.w**2
        if c2 is None:
            c2 = [0 for i in range(len(self.w))]
        if c3 is None:
            c3 = [0 for i in range(len(self.w))]

        self.id = yaeos_c.psrk(
            tc=self.tc,
            pc=self.pc,
            w=self.w,
            ngs=number_of_groups,
            g_ids=groups_ids,
            g_v=groups_ammounts,
            c1=c1,
            c2=c2,
            c3=c3,
        )
