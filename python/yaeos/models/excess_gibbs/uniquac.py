"""UNIQUAC"""

from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class UNIQUAC(GeModel):
    """UNIQUAC"""

    def __init__(self, qs, rs, aij, bij, cij, dij, eij) -> None:
        self.qs = qs
        self.rs = rs
        self.aij = aij
        self.bij = bij
        self.cij = cij
        self.dij = dij
        self.eij = eij

        self.id = yaeos_c.uniquac(qs, rs, aij, bij, cij, dij, eij)
