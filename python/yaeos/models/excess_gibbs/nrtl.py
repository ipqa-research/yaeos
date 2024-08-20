from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class NRTL(GeModel):

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        self.id = yaeos_c.nrtl(a, b, c)
