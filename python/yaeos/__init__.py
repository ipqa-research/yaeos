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
    PSRK,
    PengRobinson76,
    PengRobinson78,
    SoaveRedlichKwong,
    QMR,
    QMRTD,
    CMR,
    CMRTD,
    RKPR,
    sDDLC,
    HV,
    HVNRTL,
    MHV,
    AlphaSoave,
    AlphaRKPR,
    AlphaMathiasCopeman,
)
from yaeos.models.residual_helmholtz.multifluid import GERG2008
from yaeos.models.residual_helmholtz.saft import PCSAFT

__all__ = [
    "envelopes",
    "constants",
    "GPEC",
    "yaeos_c",
    "SoaveRedlichKwong",
    "PengRobinson76",
    "PengRobinson78",
    "RKPR",
    "PCSAFT",
    "PSRK",
    "QMR",
    "QMRTD",
    "CMR",
    "CMRTD",
    "GERG2008",
    "NRTL",
    "sDDLC",
    "UNIFACDortmund",
    "UNIFACPSRK",
    "UNIFACVLE",
    "UNIQUAC",
    "MHV",
    "HV",
    "HVNRTL",
    "AlphaSoave",
    "AlphaRKPR",
    "AlphaMathiasCopeman",
]


__version__ = importlib.metadata.version("yaeos")
