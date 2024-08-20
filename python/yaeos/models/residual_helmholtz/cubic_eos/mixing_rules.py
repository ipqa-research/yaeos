from abc import ABC, abstractmethod

from yaeos.lib import yaeos_c


class CubicMixRule(ABC):
    @abstractmethod
    def set_mixrule(self, ar_model_id):
        raise NotImplementedError


class QMR(CubicMixRule):
    def __init__(self, kij, lij):
        self.kij = kij
        self.lij = lij

    def set_mixrule(self, ar_model_id):
        yaeos_c.set_qmr(ar_model_id, self.kij, self.lij)


class MHV(CubicMixRule):
    def __init__(self, ge, q, lij=None):
        self.ge = ge
        self.q = q
        self.lij = lij

    def set_mixrule(self, ar_model_id):
        yaeos_c.set_mhv(ar_model_id, self.ge.id, self.q)
