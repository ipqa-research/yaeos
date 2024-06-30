"""CubicEoS interface
"""

from abc import ABC, abstractmethod
from functools import partial

import numpy as np
from yaeos import yaeos_c

class GeModel(ABC):
    ...

    def __del__(self):
        yaeos_c.make_available_ge_models_list(self.id)

class ArModel(ABC):

    def fugacity(self, n, v, t, dt=None, dp=None, dn=None):

        nc = len(n)

        if dt:
            dt = np.empty(nc, order="F")
        if dp:
            dp = np.empty(nc, order="F")
        if dn:
            dn = np.empty((nc, nc), order="F")

        res = yaeos_c.fug_vt(self.id, n, v, t, dlnphidt=dt, dlnphidp=dp, dlnphidn=dn)
        res = {"ln_phi": res, "dt": dt, "dp": dp, "dn": dn}
        return res

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
    def set_mixrule(ar_model_id):
        raise NotImplementedError


class CubicEoS(ArModel):
    def __init__(self, tc, pc, w):
        nc = len(tc)
        self.nc = nc
        self.tc = tc
        self.pc = pc
        self.w = w


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

    def __init__(self, tc, pc, w, mixrule: CubicMixRule):
        super(PengRobinson76, self).__init__(tc, pc, w)
        self.id = yaeos_c.pr76(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        mixrule.set_mixrule(self.id)

    def set_mixrule(self, mixrule: CubicMixRule):
        self.mixrule = mixrule
        self.mixrule.set_mixrule(self.id)

