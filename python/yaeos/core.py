"""CubicEoS interface
"""
from abc import ABC, abstractmethod
from functools import partial

import numpy as np
from yaeos import yaeos_c as yaeos


class ArModel(ABC):

    def fugacity(self, n, v, t):
        return yaeos.fug_vt(self.id, n, v, t)

    def __del__(self):
        yaeos.make_available_ar_models_list(self.id)


class PengRobinson76(ArModel):
    name = "PengRobinson76"

    def __init__(self, tc, pc, w, kij=None, lij=None):
        self.id = np.array(id(self), np.int64)

        n = len(tc)
        self.tc = tc
        self.pc = pc
        self.w = w

        if kij is None:
            self.kij = np.zeros((n, n))
        else:
            self.kij = kij

        if lij is None:
            self.lij = np.zeros((n, n))
        else:
            self.lij = lij

        self.id = yaeos.pr76(self.tc, self.pc, self.w, self.kij, self.lij)


def run_models():
    tc = np.array([190, 310])
    pc = np.array([14, 30])
    w = np.array([0.001, 0.03])

    kij = np.zeros((2, 2))
    kij[0, 1] = 0.1
    kij[1, 0] = 0.1

    lij = kij / 2

    m1 = PengRobinson76(tc, pc, w, kij, lij)
    m2 = PengRobinson76(tc / 2, pc, w, kij, lij)

    n = [0.3, 0.7]
    v = 1
    t = 150

    print(m1.fugacity(n, v, t))
    print(m2.fugacity(n, v, t))
    print(m1.fugacity(n, v, t))


# run_models()
