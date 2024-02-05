"""CubicEoS interface
"""
import numpy as np
from yaeos import yaeos_python


class CubicEoS:

    def __init__(self, tc, pc, w, kij, lij):
        self.tc = tc
        self.pc = pc
        self.w = w
        self.kij = kij
        self.lij = lij
        self.id = np.array(id(self), np.int64)

    def set_parameters(foo):
        from functools import wraps

        @wraps(foo)
        def checker(self, *args, **kwargs):
            if self.id != yaeos_python.running_model:
                yaeos_python.running_model = self.id
                yaeos_python.pr76(tc, pc, w, kij, lij)
            return foo(self, *args, **kwargs)

        return checker

    @set_parameters
    def fugacity(self, n, v, t):
        return yaeos_python.fugacity(
            n, v, t
        )


tc = np.array([190, 310])
pc = np.array([14, 30])
w = np.array([0.001, 0.03])

kij = np.zeros((2, 2))
kij[0, 1] = 0.1
kij[1, 0] = 0.1
lij = kij/2

m1 = CubicEoS(tc, pc, w, kij, lij)
m2 = CubicEoS(tc, pc, w, kij, lij)

n = [0.3, 0.7]
v = 1
t = 150

m1.fugacity(n, v, t)