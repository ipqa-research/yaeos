from abc import ABC, abstractmethod

from yaeos.core import ArModel
from yaeos.lib import yaeos_c


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
