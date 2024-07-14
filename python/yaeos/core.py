"""CubicEoS interface
"""

from abc import ABC, abstractmethod

import numpy as np

from yaeos import yaeos_c


class GeModel(ABC):
    """Excess Gibbs model."""

    def __del__(self):
        yaeos_c.make_available_ge_models_list(self.id)


class ArModel(ABC):
    """Residual Helmholtz (Ar) model"""

    def fugacity(self, n, v, t, dt=None, dp=None, dn=None):

        nc = len(n)

        if dt:
            dt = np.empty(nc, order="F")
        if dp:
            dp = np.empty(nc, order="F")
        if dn:
            dn = np.empty((nc, nc), order="F")

        res = yaeos_c.fug_vt(
            self.id, n, v, t, dlnphidt=dt, dlnphidp=dp, dlnphidn=dn
        )
        res = {"ln_phi": res, "dt": dt, "dp": dp, "dn": dn}
        return res

    def flash_tp(self, z, temperature, pressure):
        """Two-phase split with specification of temperature and pressure.

        Calculates the phase split at a given pressure and temperature
        """

        x, y, pressure, temperature, volume_x, volume_y, beta = yaeos_c.flash(
            self.id, z, temperature, pressure
        )

        flash_result = {
            "x": x,
            "y": y,
            "T": temperature,
            "P": pressure,
            "beta": beta,
        }

        return flash_result

    def saturation_pressure(self, z, temperature, kind="bubble"):
        """Saturation pressure at specified temperature

        Arguments
        ---------
        z: array
            Global molar fractions
        temperature: float
            Temperature [K]
        kind: string
            Kind of saturation point, defaults to "bubble". Options are
                - "bubble"
                - "dew"
                - "liquid-liquid"
        """
        p, x, y, volume_x, volume_y, beta = yaeos_c.saturation_pressure(
            self.id, z, temperature, kind
        )

        return {
            "x": x,
            "y": y,
            "Vx": volume_x,
            "Vy": volume_y,
            "T": temperature,
            "P": p,
            "beta": beta,
        }

    def __del__(self):
        yaeos_c.make_available_ar_models_list(self.id)


class NRTL:

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        self.id = yaeos_c.nrtl(a, b, c)


class CubicMixRule(ABC):
    @abstractmethod
    def set_mixrule(self, ar_model_id):
        raise NotImplementedError


class CubicEoS(ArModel):
    def __init__(self, tc, pc, w):
        nc = len(tc)
        self.nc = nc
        self.tc = tc
        self.pc = pc
        self.w = w

    def set_mixrule(self, mixrule: CubicMixRule):
        self.mixrule = mixrule
        self.mixrule.set_mixrule(self.id)


class QMR(ABC):

    def __init__(self, kij, lij):
        self.kij = kij
        self.lij = lij

    def set_mixrule(self, ar_model_id):
        yaeos_c.set_qmr(ar_model_id, self.kij, self.lij)


class MHV(ABC):

    def __init__(self, ge, q, lij=None):
        self.ge = ge
        self.q = q
        self.lij = lij

    def set_mixrule(self, ar_model_id):
        yaeos_c.set_mhv(ar_model_id, self.ge.id, self.q)


class PengRobinson76(CubicEoS):
    name = "PengRobinson76"

    def __init__(self, tc, pc, w, mixrule: CubicMixRule = None):
        super(PengRobinson76, self).__init__(tc, pc, w)
        self.id = yaeos_c.pr76(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)


class PengRobinson78(CubicEoS):
    name = "PengRobinson78"

    def __init__(self, tc, pc, w, mixrule: CubicMixRule = None):
        super(PengRobinson78, self).__init__(tc, pc, w)
        self.id = yaeos_c.pr78(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)


class SoaveRedlichKwong(CubicEoS):
    name = "SoaveReldichKwong"

    def __init__(self, tc, pc, w, mixrule: CubicMixRule = None):
        super(SoaveRedlichKwong, self).__init__(tc, pc, w)
        self.id = yaeos_c.srk(self.tc, self.pc, self.w)
        if mixrule:
            mixrule.set_mixrule(self.id)


class RKPR(CubicEoS):
    name = "RKPR"

    def __init__(
        self, tc, pc, w, zc, k=None, delta_1=None, mixrule: CubicMixRule = None
    ):
        super(RKPR, self).__init__(tc, pc, w)
        self.zc = zc
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
        if mixrule:
            mixrule.set_mixrule(self.id)
