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

    def lnphi_vt(self, n, v, t, dt=None, dp=None, dn=None):

        nc = len(n)

        if dt:
            dt = np.empty(nc, order="F")
        if dp:
            dp = np.empty(nc, order="F")
        if dn:
            dn = np.empty((nc, nc), order="F")

        res = yaeos_c.lnphi_vt(
            self.id, n, v, t, dlnphidt=dt, dlnphidp=dp, dlnphidn=dn
        )
        res = {"ln_phi": res, "dt": dt, "dp": dp, "dn": dn}
        return res

    def lnphi_pt(self, n, p, t, root="stable", dt=None, dp=None, dn=None):

        nc = len(n)

        if dt:
            dt = np.empty(nc, order="F")
        if dp:
            dp = np.empty(nc, order="F")
        if dn:
            dn = np.empty((nc, nc), order="F")

        res = yaeos_c.lnphi_pt(
            self.id, n, p, t, root, dlnphidt=dt, dlnphidp=dp, dlnphidn=dn
        )
        res = {"ln_phi": res, "dt": dt, "dp": dp, "dn": dn}
        return res

    def pressure(self, n, v, t, dv=None, dt=None, dn=None):
        nc = len(n)

        if dv:
            dv = np.empty(1, order="F")
        if dt:
            dt = np.empty(1, order="F")
        if dn:
            dn = np.empty(nc, order="F")

        res = yaeos_c.pressure(self.id, n, v, t, dpdv=dv, dpdt=dt, dpdn=dn)
        res = {"P": res, "dv": dv, "dt": dt, "dn": dn}
        return res

    def volume(self, n, p, t, root="stable"):

        res = yaeos_c.volume(self.id, n, p, t, root)
        res = {"V": res}
        return res

    def flash_pt(self, z, pressure, temperature):
        """Two-phase split with specification of temperature and pressure.

        Calculates the phase split at a given pressure and temperature
        """

        x, y, pressure, temperature, volume_x, volume_y, beta = yaeos_c.flash(
            self.id, z, p=pressure, t=temperature
        )

        flash_result = {
            "x": x,
            "y": y,
            "P": pressure,
            "T": temperature,
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

    def phase_envelope_pt(
        self, z, kind="bubble", max_points=300, T0=150, P0=1
    ):
        Ts, Ps, Tcs, Pcs = yaeos_c.pt2_phase_envelope(
            self.id, z, kind=kind, t0=T0, p0=P0, max_points=max_points
        )
        return Ts, Ps, Tcs, Pcs

    def __del__(self):
        yaeos_c.make_available_ar_models_list(self.id)


class NRTL:

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        self.id = yaeos_c.nrtl(a, b, c)
