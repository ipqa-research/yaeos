"""Yet Another Equation-Of-State (library).

Library to use EoS-based calculations. This main module imports all the
relevant constants, procedures and objects to have better access to them.
"""

__version__ = "1.5.0"

from yaeos.lib import yaeos_c
from yaeos.models.excess_gibbs import NRTL, UNIFACVLE, UNIQUAC
from yaeos.models.residual_helmholtz.cubic_eos import (
    HV,
    MHV,
    PSRK,
    PengRobinson76,
    PengRobinson78,
    QMR,
    QMRTD,
    RKPR,
    SoaveRedlichKwong,
)


__all__ = [
    "yaeos_c",
    "SoaveRedlichKwong",
    "PengRobinson76",
    "PengRobinson78",
    "RKPR",
    "PSRK",
    "QMR",
    "QMRTD",
    "NRTL",
    "UNIFACVLE",
    "UNIQUAC",
    "MHV",
    "HV",
]
