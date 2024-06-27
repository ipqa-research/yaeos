"""CubicEoS interface
"""
from abc import ABC, abstractmethod
from functools import partial

import numpy as np
from yaeos import yaeos_c as yaeos


class ArModel(ABC):

    def fugacity(self, n, v, t, dt=None, dp=None, dn=None):

        nc = len(n)

        if dt:
            dt = np.empty(nc, order="F")
        if dp:
            dp = np.empty(nc, order="F")
        if dn:
            dn = np.empty((nc, nc), order="F")

        res = yaeos.fug_vt(
            self.id, 
            n, v, t, 
            dlnphidt=dt, dlnphidp=dp, dlnphidn=dn
            )
        res = {
            "ln_phi": res, "dt": dt, "dp": dp, "dn": dn
        }
        return res

    def __del__(self):
        yaeos.make_available_ar_models_list(self.id)


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
        yaeos.set_qmr(ar_model_id, self.kij, self.lij)


class PengRobinson76(CubicEoS):
    name = "PengRobinson76"
    
    def __init__(self, tc, pc, w, mixrule: CubicMixRule):
        super(PengRobinson76, self).__init__(tc, pc, w)
        self.id = yaeos.pr76(self.tc, self.pc, self.w)
        mixrule.set_mixrule(self.id)