"""Cubic EoS module.

Implemented models:

- Cubic EoS:
    - PengRobinson76: Peng-Robinson model (1976)
    - PengRobinson78: Peng-Robinson model (1978)
    - SoaveRedlichKwong: Soave-Redlich-Kwong model
    - RKPR: RKPR model
- Mixing rules: mixing rules for cubic EoS
    - QMR: cuadratic mixing rule
    - MHV: modified Huron-Vidal mixing rule
"""

from .cubic_eos import (
    CubicEoS,
    PSRK,
    PengRobinson76,
    PengRobinson78,
    RKPR,
    SoaveRedlichKwong,
)
from .mixing_rules import CubicMixRule, HV, HVNRTL, MHV, QMR, QMRTD


__all__ = [
    "CubicEoS",
    "PengRobinson76",
    "PengRobinson78",
    "SoaveRedlichKwong",
    "RKPR",
    "PSRK",
    "CubicMixRule",
    "QMR",
    "QMRTD",
    "MHV",
    "HV",
    "HVNRTL",
]
