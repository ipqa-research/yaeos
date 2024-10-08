"""Residual Helmholtz energy models module.

Yaeos Residual Helmholtz module. This module provides the following submodules:

- residual_helmholtz: Residual Helmholtz energy models
    - Cubic EoS:
        - PengRobinson76: Peng-Robinson model (1976)
        - PengRobinson78: Peng-Robinson model (1978)
        - SoaveRedlichKwong: Soave-Redlich-Kwong model
        - RKPR: RKPR model
        - PSRK: Predictive-Soave-Redlich-Kwong model
    - Mixing rules: mixing rules for cubic EoS
        - QMR: cuadratic mixing rule
        - MHV: modified Huron-Vidal mixing rule
"""

from . import cubic_eos


__all__ = ["cubic_eos"]
