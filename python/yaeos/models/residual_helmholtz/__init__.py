"""Residual Helmholtz energy models module.

Yaeos Residual Helmholtz module. This module provides the following submodules:

- residual_helmholtz: Residual Helmholtz energy models
    - Cubic EoS:
        - Cubic EoS:
            - PengRobinson76: Peng-Robinson model (1976)
            - PengRobinson78: Peng-Robinson model (1978)
            - SoaveRedlichKwong: Soave-Redlich-Kwong model
            - RKPR: RKPR model
            - PSRK: Predictive-Soave-Redlich-Kwong model
        - Mixing rules: mixing rules for cubic EoS
            - QMR: cuadratic mixing rule
            - HV: Huron-Vidal mixing rule
            - MHV: Michelsen's modified Huron-Vidal mixing rule
            - HVNRTL: Huron-Vidal NRTL mixing rule
    - Multifluids Equations of State (Multifluid EoS):
        - GERG2008: GERG2008 Residual contribution
    - SAFT Equations of State (SAFT EoS):
        - PC-SAFT: PC-SAFT Residual contribution
"""

from yaeos.models.residual_helmholtz import cubic_eos
from yaeos.models.residual_helmholtz import multifluid
from yaeos.models.residual_helmholtz import saft

__all__ = ["cubic_eos", "multifluid", "saft"]
