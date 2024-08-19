"""yaeos Python API core module.

ArModel and GeModel abstract classes definition. Also, the implementation of 
the models' thermoprops methods.
"""

from abc import ABC
from typing import Union

import numpy as np

from yaeos.lib import yaeos_c


class GeModel(ABC):
    """Excess Gibbs (Ge) model abstract class."""

    def __del__(self):
        yaeos_c.make_available_ge_models_list(self.id)


class ArModel(ABC):
    """Residual Helmholtz (Ar) model abstract class."""

    def lnphi_vt(
        self,
        moles,
        volume: float,
        temperature: float,
        dt: bool = False,
        dp: bool = False,
        dn: bool = False,
    ) -> Union[np.ndarray, tuple[np.ndarray, dict]]:
        """Calculate fugacity coefficent given volume and temperature.

        Calculate :math:`ln \phi_i(n,V,T)` and its derivatives with respect to
        temperature, pressure and moles number.

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]
        dt : bool, optional
            Calculate temperature derivative, by default False
        dp : bool, optional
            Calculate pressure derivative, by default False
        dn : bool, optional
            Calculate moles derivative, by default False

        Returns
        -------
        Union[np.ndarray, tuple[np.ndarray, dict]]
            :math:`ln \phi_i(n,V,T)` vector or tuple with
            :math:`ln \phi_i(n,V,T)` vector and derivatives dictionary if any
            derivative is asked

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating ln_phi only
            # will print: [-1.45216274 -2.01044828]

            print(model.lnphi_vt([5.0, 5.6], 1.0, 300.0))

            # Asking for derivatives
            # will print:
            # (
            # array([-1.45216274, -2.01044828]),
            # {'dt': array([0.01400063, 0.01923493]), 'dp': None, 'dn': None}
            # )

            print(model.lnphi_vt([5.0, 5.6], 1.0, 300.0, dt=True)
        """
        nc = len(moles)

        dt = np.empty(nc, order="F") if dt else None
        dp = np.empty(nc, order="F") if dp else None
        dn = np.empty((nc, nc), order="F") if dn else None

        res = yaeos_c.lnphi_vt(
            self.id,
            moles,
            volume,
            temperature,
            dlnphidt=dt,
            dlnphidp=dp,
            dlnphidn=dn,
        )

        if dt is None and dp is None and dn is None:
            ...
        else:
            res = (res, {"dt": dt, "dp": dp, "dn": dn})
        return res

    def lnphi_pt(
        self,
        moles,
        pressure: float,
        temperature: float,
        root: str = "stable",
        dt: bool = False,
        dp: bool = False,
        dn: bool = False,
    ) -> Union[np.ndarray, tuple[np.ndarray, dict]]:
        """Calculate fugacity coefficent given pressure and temperature.

        Calculate :math:`ln \phi_i(n,P,T)` and its derivatives with respect to
        temperature, pressure and moles number.

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        pressure : float
            Pressure [bar]
        temperature : float
            Temperature [K]
        root : str, optional
            Volume root, use: "liquid", "vapor" or "stable", by default
            "stable"
        dt : bool, optional
            Calculate temperature derivative, by default False
        dp : bool, optional
            Calculate pressure derivative, by default False
        dn : bool, optional
            Calculate moles derivative, by default False

        Returns
        -------
        Union[np.ndarray, tuple[np.ndarray, dict]]
            :math:`ln \phi_i(n,P,T)` vector or tuple with
            :math:`ln \phi_i(n,P,T)` vector and derivatives dictionary if any
            derivative is asked

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating ln_phi only
            # will print: [-0.10288733 -0.11909807]

            print(model.lnphi_pt([5.0, 5.6], 10.0, 300.0))

            # Asking for derivatives
            # will print:
            # (
            # array([-0.10288733, -0.11909807]),
            # {'dt': array([0.00094892, 0.00108809]), 'dp': None, 'dn': None}
            # )

            print(model.lnphi_pt([5.0, 5.6], 10.0, 300.0, dt=True)
        """
        nc = len(moles)

        dt = np.empty(nc, order="F") if dt else None
        dp = np.empty(nc, order="F") if dp else None
        dn = np.empty((nc, nc), order="F") if dn else None

        res = yaeos_c.lnphi_pt(
            self.id,
            moles,
            pressure,
            temperature,
            root,
            dlnphidt=dt,
            dlnphidp=dp,
            dlnphidn=dn,
        )

        if dt is None and dp is None and dn is None:
            ...
        else:
            res = (res, {"dt": dt, "dp": dp, "dn": dn})
        return res

    def pressure(
        self,
        moles,
        volume: float,
        temperature: float,
        dv: bool = False,
        dt: bool = False,
        dn: bool = False,
    ) -> Union[float, tuple[float, dict]]:
        """Calculate pressure given volume and temperature [bar].

        Calculate :math:`P(n,V,T)` and its derivatives with respect to
        volume, temperature and moles number.

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]
        dv : bool, optional
            Calculate volume derivative, by default False
        dt : bool, optional
            Calculate temperature derivative, by default False
        dn : bool, optional
            Calculate moles derivative, by default False

        Returns
        -------
        Union[float, tuple[float, dict]]
            :math:`P(n,V,T)` or tuple with :math:`P(n,V,T)` and derivatives
            dictionary if any derivative is asked [bar]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating pressure only
            # will print: 16.011985733846956

            print(model.pressure(np.array([5.0, 5.6]), 2.0, 300.0))

            # Asking for derivatives
            # will print:
            # (
            # 16.011985733846956,
            # {'dv': None, 'dt': np.float64(0.7664672352866752), 'dn': None}
            # )

            print(model.pressure(np.array([5.0, 5.6]), 2.0, 300.0, dt=True))
        """
        nc = len(moles)

        dv = np.empty(1, order="F") if dv else None
        dt = np.empty(1, order="F") if dt else None
        dn = np.empty(nc, order="F") if dn else None

        res = yaeos_c.pressure(
            self.id, moles, volume, temperature, dpdv=dv, dpdt=dt, dpdn=dn
        )

        if dt is None and dv is None and dn is None:
            ...
        else:
            res = (
                res,
                {
                    "dv": dv if dv is None else dv[0],
                    "dt": dt if dt is None else dt[0],
                    "dn": dn,
                },
            )
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
            "Vx": volume_x,
            "Vy": volume_y,
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
        self, z, kind="bubble", max_points=300, t0=150, p0=1
    ):
        ts, ps, tcs, pcs = yaeos_c.pt2_phase_envelope(
            self.id, z, kind=kind, t0=t0, p0=p0, max_points=max_points
        )
        return ts, ps, tcs, pcs

    def __del__(self):
        yaeos_c.make_available_ar_models_list(self.id)


class NRTL:

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        self.id = yaeos_c.nrtl(a, b, c)
