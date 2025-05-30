"""Yet Another Equation-Of-State (library).

Library to use EoS-based calculations. This main module imports all the
relevant constants, procedures and objects to have better access to them.
"""

import importlib.metadata

import yaeos.constants as constants
import yaeos.envelopes as envelopes
from yaeos.gpec import GPEC
from yaeos.lib import yaeos_c
from yaeos.models.excess_gibbs import (
    NRTL,
    UNIFACDortmund,
    UNIFACPSRK,
    UNIFACVLE,
    UNIQUAC,
)
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
    "envelopes",
    "constants",
    "GPEC",
    "yaeos_c",
    "SoaveRedlichKwong",
    "PengRobinson76",
    "PengRobinson78",
    "RKPR",
    "PSRK",
    "QMR",
    "QMRTD",
    "NRTL",
    "UNIFACDortmund",
    "UNIFACPSRK",
    "UNIFACVLE",
    "UNIQUAC",
    "MHV",
    "HV",
]


__version__ = importlib.metadata.version("yaeos")
