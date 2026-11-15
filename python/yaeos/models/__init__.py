"""Models module.

Yaeos models module. This module provides the following submodules:

- excess_gibbs: Excess Gibbs energy models
    - NRTL: non-random two-liquid model
    - UNIFACVLE: Original UNIFAC VLE model
    - UNIFACPSRK: PSRK-UNIFAC model
    - UNIFACDortmund: UNIFAC Dortmund model
    - UNIQUAC: UNIversal QUAsiChemical Excess Gibbs free energy model

- residual_helmholtz: Residual Helmholtz energy models
    - Cubic EoS:
        - PengRobinson76: Peng-Robinson model (1976)
        - PengRobinson78: Peng-Robinson model (1978)
        - SoaveRedlichKwong: Soave-Redlich-Kwong model
        - RKPR: RKPR model
        - PSRK: Predictive Soave-Redlich-Kwong model
    - Alpha functions: AlphaSoave, AlphaRKPR, AlphaMathiasCopeman
    - Mixing rules: mixing rules for cubic EoS
        - QMR: quadratic mixing rule
        - QMRTD: quadratic mixing rule with temperature dependence
        - CMR: cubic mixing rule
        - CMRTD: cubic mixing rule with temperature dependence
        - sDDLC: segmented density dependent local composition mixing rule
        - MHV: modified Huron-Vidal mixing rule
        - HV: Huron-Vidal mixing rule
        - HVNRTL: Huron-Vidal NRTL mixing rule
    - Multifluid EoS:
        - GERG2008: GERG2008 Residual contribution
    - SAFT EoS:
        - PC-SAFT: PC-SAFT Residual contribution
"""

from yaeos.models import excess_gibbs, residual_helmholtz

__all__ = ["excess_gibbs", "residual_helmholtz"]
