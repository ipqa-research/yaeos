"""Models module.

Yaeos models module. This module provides the following submodules:

- excess_gibbs: Excess Gibbs energy models
    - NRTL: non-random two-liquid model
    - UNIFACVLE: Original UNIFAC VLE model
    - UNIQUAC: UNIversal QUAsiChemical Excess Gibbs free energy model

- residual_helmholtz: Residual Helmholtz energy models
    - Cubic EoS:
        - PengRobinson76: Peng-Robinson model (1976)
        - PengRobinson78: Peng-Robinson model (1978)
        - SoaveRedlichKwong: Soave-Redlich-Kwong model
        - RKPR: RKPR model
    - Mixing rules: mixing rules for cubic EoS
        - QMR: cuadratic mixing rule
        - MHV: modified Huron-Vidal mixing rule
"""

from . import excess_gibbs, residual_helmholtz


__all__ = ["excess_gibbs", "residual_helmholtz"]
