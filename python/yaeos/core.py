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

    def ln_gamma(
        self, moles, temperature: float, dt: bool = False, dn: bool = False
    ) -> Union[np.ndarray, tuple[np.ndarray, dict]]:
        r"""Calculate natural logarithm of activity coefficients.

        Calculate :math:`\ln \gamma_i(n,T)` vector.

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        temperature : float
            Temperature [K]

        Returns
        -------
        np.ndarray
            :math:`ln \gamma_i(n,T)` vector

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import NRTL


            a = np.array([[0, 0.3], [0.3, 0]])
            b = np.array([[0, 0.4], [0.4, 0]])
            c = np.array([[0, 0.5], [0.5, 0]])

            nrtl = NRTL(a, b, c)

            # Evaluating ln_gamma only
            print(nrtl.ln_gamma([5.0, 5.6], 300.0))

            # Asking for derivatives

            print(nrtl.ln_gamma([5.0, 5.6], 300.0, dt=True, dn=True))
        """
        nc = len(moles)

        dt = np.empty(nc, order="F") if dt else None
        dn = np.empty((nc, nc), order="F") if dn else None

        res = yaeos_c.ln_gamma_ge(
            self.id,
            moles,
            temperature,
            dlngamma_dt=dt,
            dlngamma_dn=dn,
        )

        if dt is None and dn is None:
            ...
        else:
            res = (res, {"dt": dt, "dn": dn})
        return res

    def excess_gibbs(
        self,
        moles,
        temperature: float,
        dt: bool = False,
        dt2: bool = False,
        dn: bool = False,
        dtn: bool = False,
        dn2: bool = False,
    ) -> Union[np.ndarray, tuple[np.ndarray, dict]]:
        """Calculate excess Gibbs energy [bar L].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        temperature : float
            Temperature [K]
        dt : bool, optional
            Calculate temperature derivative, by default False
        dt2 : bool, optional
            Calculate temperature second derivative, by default False
        dn : bool, optional
            Calculate moles derivative, by default False
        dtn : bool, optional
            Calculate cross temperature and moles derivative, by default False
        dn2 : bool, optional
            Calculate moles second derivative, by default False

        Returns
        -------
        Union[np.ndarray, tuple[np.ndarray, dict]]
            Excess Gibbs energy or tuple with excess Gibbs energy and
            derivatives dictionary if any derivative is asked [bar L]

        Example
        -------
        .. code-block:: python

            from yaeos import UNIFACVLE

            # Ethanol - water system
            groups = [{1: 2, 2: 1, 14: 1}, {16: 1}]

            model = UNIFACVLE(groups)

            # Evaluating excess Gibbs energy only
            print(model.excess_gibbs(model.excess_gibbs([0.5,0.5], 303.15))

            # Asking for derivatives
            print(
                model.excess_gibbs(
                    [0.5,0.5],
                    303.15,
                    dt=True,
                    dt2=True,
                    dn=True,
                    dtn=True,
                    dn2=True
            )
        """
        nc = len(moles)

        dt = np.empty(1, order="F") if dt else None
        dt2 = np.empty(1, order="F") if dt2 else None
        dn = np.empty(nc, order="F") if dn else None
        dtn = np.empty(nc, order="F") if dtn else None
        dn2 = np.empty((nc, nc), order="F") if dn2 else None

        all_none = (
            dt is None
            and dt2 is None
            and dn is None
            and dtn is None
            and dn2 is None
        )

        res = yaeos_c.excess_gibbs_ge(
            self.id,
            moles,
            temperature,
            get=dt,
            get2=dt2,
            gen=dn,
            getn=dtn,
            gen2=dn2,
        )

        if all_none:
            ...
        else:
            res = (
                res,
                {
                    "dt": dt if dt is None else dt[0],
                    "dt2": dt2 if dt2 is None else dt2[0],
                    "dn": dn,
                    "dtn": dtn,
                    "dn2": dn2,
                },
            )

        return res

    def excess_enthalpy(
        self, moles, temperature: float, dt: bool = False, dn: bool = False
    ) -> Union[np.ndarray, tuple[np.ndarray, dict]]:
        """Calculate excess enthalpy [bar L].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        temperature : float
            Temperature [K]
        dt : bool, optional
            Calculate temperature derivative, by default False
        dn : bool, optional
            Calculate moles derivative, by default False

        Returns
        -------
        Union[np.ndarray, tuple[np.ndarray, dict]]
            Excess enthalpy or tuple with excess enthalpy and derivatives
            dictionary if any derivative is asked [bar L]

        Example
        -------
        .. code-block:: python

            from yaeos import UNIFACVLE

            # Ethanol - water system
            groups = [{1: 2, 2: 1, 14: 1}, {16: 1}]

            model = UNIFACVLE(groups)

            # Evaluating excess enthalpy only
            print(model.excess_enthalpy([0.5, 0.5], 303.15))

            # Asking for derivatives
            print(model.excess_enthalpy([0.5, 0.5], 303.15, dt=True, dn=True))
        """
        nc = len(moles)

        dt = np.empty(1, order="F") if dt else None
        dn = np.empty(nc, order="F") if dn else None

        res = yaeos_c.excess_enthalpy_ge(
            self.id,
            moles,
            temperature,
            het=dt,
            hen=dn,
        )

        if dt is None and dn is None:
            ...
        else:
            res = (res, {"dt": dt if dt is None else dt[0], "dn": dn})

        return res

    def excess_entropy(
        self, moles, temperature: float, dt: bool = False, dn: bool = False
    ) -> Union[np.ndarray, tuple[np.ndarray, dict]]:
        """Calculate excess entropy [bar L / K].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        temperature : float
            Temperature [K]
        dt : bool, optional
            Calculate temperature derivative, by default False
        dn : bool, optional
            Calculate moles derivative, by default False

        Returns
        -------
        Union[np.ndarray, tuple[np.ndarray, dict]]
            Excess entropy or tuple with excess entropy and derivatives
            dictionary if any derivative is asked [bar L / K]

        Example
        -------
        .. code-block:: python

            from yaeos import UNIFACVLE

            # Ethanol - water system
            groups = [{1: 2, 2: 1, 14: 1}, {16: 1}]

            model = UNIFACVLE(groups)

            # Evaluating excess entropy only
            print(model.excess_entropy([0.5, 0.5], 303.15))

            # Asking for derivatives
            print(model.excess_entropy([0.5, 0.5], 303.15, dt=True, dn=True))
        """
        nc = len(moles)

        dt = np.empty(1, order="F") if dt else None
        dn = np.empty(nc, order="F") if dn else None

        res = yaeos_c.excess_entropy_ge(
            self.id,
            moles,
            temperature,
            set=dt,
            sen=dn,
        )

        if dt is None and dn is None:
            ...
        else:
            res = (res, {"dt": dt if dt is None else dt[0], "dn": dn})

        return res

    def __del__(self) -> None:
        """Delete the model from the available models list (Fortran side)."""
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
        r"""Calculate fugacity coefficent given volume and temperature.

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
        r"""Calculate fugacity coefficent given pressure and temperature.

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
            Pressure or tuple with Presure and derivatives dictionary if any
            derivative is asked [bar]

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

    def volume(
        self, moles, pressure: float, temperature: float, root: str = "stable"
    ) -> float:
        """Calculate volume given pressure and temperature [L].

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

        Returns
        -------
        float
            Volume [L]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating stable root volume
            # will print: 23.373902973572587

            print(model.volume(np.array([5.0, 5.6]), 10.0, 300.0))

            # Liquid root volume (not stable)
            # will print: 0.8156388756398074

            print(model.volume(np.array([5.0, 5.6]), 10.0, 300.0, "liquid"))
        """
        res = yaeos_c.volume(self.id, moles, pressure, temperature, root)
        return res

    def enthalpy_residual_vt(
        self,
        moles,
        volume: float,
        temperature: float,
        dt: bool = False,
        dv: bool = False,
        dn: bool = False,
    ) -> Union[float, tuple[float, dict]]:
        """Calculate residual enthalpy given volume and temperature [bar L].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]

        Returns
        -------
        Union[float, tuple[float, dict]]
            Residual enthalpy or tuple with Residual enthalpy and derivatives
            dictionary if any derivative is asked [bar L]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating residual enthalpy only
            # will print: -182.50424367123696

            print(
                model.enthalpy_residual_vt(np.array([5.0, 5.6]), 10.0, 300.0)
            )

            # Asking for derivatives
            # will print:
            # (
            # -182.50424367123696,
            # {'dt': 0.21542452742588686, 'dv': None, 'dn': None}
            # )

            print(
                model.enthalpy_residual_vt(
                    np.array([5.0, 5.6]),
                    10.0,
                    300.0,
                    dt=True)
                )
            )
        """
        nc = len(moles)

        dt = np.empty(1, order="F") if dt else None
        dv = np.empty(1, order="F") if dv else None
        dn = np.empty(nc, order="F") if dn else None

        res = yaeos_c.enthalpy_residual_vt(
            self.id,
            moles,
            volume,
            temperature,
            hrt=dt,
            hrv=dv,
            hrn=dn,
        )

        if dt is None and dv is None and dn is None:
            ...
        else:
            res = (
                res,
                {
                    "dt": dt if dt is None else dt[0],
                    "dv": dv if dv is None else dv[0],
                    "dn": dn,
                },
            )
        return res

    def gibbs_residual_vt(
        self,
        moles,
        volume: float,
        temperature: float,
        dt: bool = False,
        dv: bool = False,
        dn: bool = False,
    ) -> Union[float, tuple[float, dict]]:
        """Calculate residual Gibbs energy at volume and temperature [bar L].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]

        Returns
        -------
        Union[float, tuple[float, dict]]
            Residual Gibbs energy or tuple with Residual Gibbs energy and
            derivatives dictionary if any derivative is asked [bar L]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating residual gibbs energy only
            # will print: -138.60374582274

            print(model.gibbs_residual_vt(np.array([5.0, 5.6]), 10.0, 300.0))

            # Asking for derivatives
            # will print:
            # (
            # -138.60374582274,
            # {'dt': 0.289312908265414, 'dv': None, 'dn': None}
            # )

            print(
                model.gibbs_residual_vt(
                    np.array([5.0, 5.6]),
                    10.0,
                    300.0,
                    dt=True
                )
            )
        """
        nc = len(moles)

        dt = np.empty(1, order="F") if dt else None
        dv = np.empty(1, order="F") if dv else None
        dn = np.empty(nc, order="F") if dn else None

        res = yaeos_c.gibbs_residual_vt(
            self.id,
            moles,
            volume,
            temperature,
            grt=dt,
            grv=dv,
            grn=dn,
        )

        if dt is None and dv is None and dn is None:
            ...
        else:
            res = (
                res,
                {
                    "dt": dt if dt is None else dt[0],
                    "dv": dv if dv is None else dv[0],
                    "dn": dn,
                },
            )
        return res

    def entropy_residual_vt(
        self,
        moles,
        volume: float,
        temperature: float,
        dt: bool = False,
        dv: bool = False,
        dn: bool = False,
    ) -> Union[float, tuple[float, dict]]:
        """Calculate residual entropy given volume and temperature [bar L / K].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]

        Returns
        -------
        Union[float, tuple[float, dict]]
            Residual entropy or tuple with Residual entropy and derivatives
            dictionary if any derivative is asked [bar L]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating residual entropy only
            # will print: -0.1463349928283233

            print(model.entropy_residual_vt(np.array([5.0, 5.6]), 10.0, 300.0))

            # Asking for derivatives
            # will print:
            # (
            # (-0.1463349928283233,
            # {'dt': 0.00024148870662932045, 'dv': None, 'dn': None})
            # )

            print(
                model.entropy_residual_vt(
                    np.array([5.0, 5.6]),
                    10.0,
                    300.0,
                    dt=True
                )
            )
        """
        nc = len(moles)

        dt = np.empty(1, order="F") if dt else None
        dv = np.empty(1, order="F") if dv else None
        dn = np.empty(nc, order="F") if dn else None

        res = yaeos_c.entropy_residual_vt(
            self.id,
            moles,
            volume,
            temperature,
            srt=dt,
            srv=dv,
            srn=dn,
        )

        if dt is None and dv is None and dn is None:
            ...
        else:
            res = (
                res,
                {
                    "dt": dt if dt is None else dt[0],
                    "dv": dv if dv is None else dv[0],
                    "dn": dn,
                },
            )
        return res

    def cv_residual_vt(
        self, moles, volume: float, temperature: float
    ) -> float:
        """Residual isochoric heat capacity given V and T [bar L / K].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]

        Returns
        -------
        float
            Residual isochoric heat capacity [bar L / K]

        Example
        -------
            .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating residual isochoric heat capacity only
            # will print: 0.07244661198879614

            print(model.cv_residual_vt(np.array([5.0, 5.6]), 10.0, 300.0))
        """
        return yaeos_c.cv_residual_vt(self.id, moles, volume, temperature)

    def cp_residual_vt(
        self, moles, volume: float, temperature: float
    ) -> float:
        """Calculate residual isobaric heat capacity given V and T [bar L / K].

        Parameters
        ----------
        moles : array_like
            Moles number vector [mol]
        volume : float
            Volume [L]
        temperature : float
            Temperature [K]

        Returns
        -------
        float
            Residual isochoric heat capacity [bar L / K]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([320.0, 375.0])   # critical temperatures [K]
            pc = np.array([45.0, 60.0])     # critical pressures [bar]
            w = np.array([0.0123, 0.045])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Evaluating residual isobaric heat capacity only
            # will print: 1.4964025088916886

            print(model.cp_residual_vt(np.array([5.0, 5.6]), 10.0, 300.0))
        """
        return yaeos_c.cp_residual_vt(self.id, moles, volume, temperature)

    def pure_saturation_pressures(
        self, component, stop_pressure=0.01, stop_temperature=100
    ):
        """Calculate pure component saturation pressures [bar].

        Calculation starts from the critical point and goes down to the
        stop pressure or stop temperature.

        Parameters
        ----------
        component : int
            Component index (starting from 1)
        stop_pressure : float, optional
            Stop pressure [bar], by default 0.01
        stop_temperature : float, optional
            Stop temperature [K], by default 100

        Returns
        -------
        dict
            Pure component saturation points dictionary with keys:
                - T: Temperature [K]
                - P: Pressure [bar]
                - Vx: Liquid Phase Volume [L/mole]
                - Vy: Vapor Phase Volume [L/mole]

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76

            tc = np.array([320.0, 375.0])
            pc = np.array([45.0, 60.0])
            w = np.array([0.0123, 0.045])

            model = PengRobinson76(tc, pc, w)
        """
        p, t, vx, vy = yaeos_c.pure_saturation_line(
            self.id, component, stop_pressure, stop_temperature
        )

        msk = ~np.isnan(t)

        return {"T": t[msk], "P": p[msk], "Vx": vx[msk], "Vy": vy[msk]}

    def flash_pt(
        self, z, pressure: float, temperature: float, k0=None
    ) -> dict:
        """Two-phase split with specification of temperature and pressure.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        pressure : float
            Pressure [bar]
        temperature : float
            Temperature [K]
        k0 : array_like, optional
            Initial guess for the split, by default None (will use k_wilson)

        Returns
        -------
        dict
            Flash result dictionary with keys:
                - x: heavy phase mole fractions
                - y: light phase mole fractions
                - Vx: heavy phase volume [L]
                - Vy: light phase volume [L]
                - P: pressure [bar]
                - T: temperature [K]
                - beta: light phase fraction

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([369.83, 507.6])       # critical temperatures [K]
            pc = np.array([42.48, 30.25])        # critical pressures [bar]
            w = np.array([0.152291, 0.301261])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Flash calculation
            # will print:
            # {
            #   'x': array([0.3008742, 0.6991258]),
            #   'y': array([0.85437317, 0.14562683]),
            #   'Vx': 0.12742569165483714,
            #   'Vy': 3.218831515959867,
            #   'P': 8.0,
            #   'T': 350.0,
            #   'beta': 0.35975821044266726
            # }

            print(model.flash_pt([0.5, 0.5], 8.0, 350.0))
        """
        if k0 is None:
            k0 = [0 for i in range(len(z))]
        x, y, pressure, temperature, volume_x, volume_y, beta = yaeos_c.flash(
            self.id, z, p=pressure, t=temperature, k0=k0
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

    def flash_pt_grid(
        self, z, pressures, temperatures, parallel=False
    ) -> dict:
        """Two-phase split with specification of temperature and pressure grid.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        pressures : array_like
            Pressures grid [bar]
        temperatures : array_like
            Temperatures grid [K]
        parallel : bool, optional
            Use parallel processing, by default False

        Returns
        -------
        dict
            Flash grid result dictionary with keys:
                - x: heavy phase mole fractions
                - y: light phase mole fractions
                - Vx: heavy phase volume [L]
                - Vy: light phase volume [L]
                - P: pressure [bar]
                - T: temperature [K]
                - beta: light phase fraction

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76

            tc = np.array([369.83, 507.6])       # critical temperatures [K]
            pc = np.array([42.48, 30.25])        # critical pressures [bar]
            w = np.array([0.152291, 0.301261])

            temperatures = [350.0, 360.0, 400.0]
            pressures = [10, 20, 30]
        """
        xs, ys, vxs, vys, betas = yaeos_c.flash_grid(
            self.id, z, pressures, temperatures, parallel=parallel
        )

        flash = {
            "x": xs,
            "y": ys,
            "Vx": vxs,
            "Vy": vys,
            "P": pressures,
            "T": temperatures,
            "beta": betas,
        }

        return flash

    def saturation_pressure(
        self, z, temperature: float, kind: str = "bubble", p0: float = 0
    ) -> dict:
        """Saturation pressure at specified temperature.

        Arguments
        ---------
        z: array_like
            Global molar fractions
        temperature: float
            Temperature [K]
        kind: str, optional
            Kind of saturation point, defaults to "bubble". Options are
                - "bubble"
                - "dew"
                - "liquid-liquid"

        Returns
        -------
        dict
            Saturation pressure calculation result dictionary with keys:
                - x: heavy phase mole fractions
                - y: light phase mole fractions
                - Vx: heavy phase volume [L]
                - Vy: light phase volume [L]
                - P: pressure [bar]
                - T: temperature [K]
                - beta: light phase fraction

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([369.83, 507.6])       # critical temperatures [K]
            pc = np.array([42.48, 30.25])        # critical pressures [bar]
            w = np.array([0.152291, 0.301261])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Saturation pressure calculation
            # will print:
            # {
            # 'x': array([0.5, 0.5]),
            # 'y': array([0.9210035 , 0.07899651]),
            # 'Vx': 0.11974125553488875,
            # 'Vy': 1.849650524323853,
            # 'T': 350.0,
            # 'P': 12.990142036059941,
            # 'beta': 0.0
            # }

            print(model.saturation_pressure(np.array([0.5, 0.5]), 350.0))
        """
        p, x, y, volume_x, volume_y, beta = yaeos_c.saturation_pressure(
            id=self.id, z=z, t=temperature, kind=kind, p0=p0
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

    def saturation_temperature(
        self, z, pressure: float, kind: str = "bubble", t0: float = 0
    ) -> dict:
        """Saturation temperature at specified pressure.

        Arguments
        ---------
        z: array_like
            Global molar fractions
        pressure: float
            Pressure [bar]
        kind: str, optional
            Kind of saturation point, defaults to "bubble". Options are
                - "bubble"
                - "dew"
                - "liquid-liquid"

        Returns
        -------
        dict
            Saturation temperature calculation result dictionary with keys:
                - x: heavy phase mole fractions
                - y: light phase mole fractions
                - Vx: heavy phase volume [L]
                - Vy: light phase volume [L]
                - P: pressure [bar]
                - T: temperature [K]
                - beta: light phase fraction

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([369.83, 507.6])       # critical temperatures [K]
            pc = np.array([42.48, 30.25])        # critical pressures [bar]
            w = np.array([0.152291, 0.301261])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Saturation pressure calculation
            # will print:
            # {
            # 'x': array([0.5, 0.5]),
            # 'y': array([0.9210035 , 0.07899651]),
            # 'Vx': 0.11974125553488875,
            # 'Vy': 1.849650524323853,
            # 'T': 350.0,
            # 'P': 12.99,
            # 'beta': 0.0
            # }

            print(model.saturation_temperature(np.array([0.5, 0.5]), 12.99))
        """
        t, x, y, volume_x, volume_y, beta = yaeos_c.saturation_temperature(
            id=self.id, z=z, p=pressure, kind=kind, t0=t0
        )

        return {
            "x": x,
            "y": y,
            "Vx": volume_x,
            "Vy": volume_y,
            "T": t,
            "P": pressure,
            "beta": beta,
        }

    # ==============================================================
    # Phase envelopes
    # --------------------------------------------------------------
    def phase_envelope_pt(
        self,
        z,
        kind: str = "bubble",
        max_points: int = 300,
        t0: float = 150.0,
        p0: float = 1.0,
    ) -> dict:
        """Two phase envelope calculation (PT).

        Parameters
        ----------
        z : array_like
            Global mole fractions
        kind : str, optional
            Kind of saturation point to start the envelope calculation,
            defaults to "bubble". Options are
            - "bubble"
            - "dew"
        max_points : int, optional
            Envelope's maximum points to calculate (T, P), by default 300
        t0 : float, optional
            Initial guess for temperature [K] for the saturation point of kind:
            `kind`, by default 150.0
        p0 : float, optional
            Initial guess for pressure [bar] for the saturation point of kind:
            `kind`, by default 1.0

        Returns
        -------
        dict
            Envelope calculation result dictionary with keys:
                - Ts: temperatures [K]
                - Ps: pressures [bar]
                - Tcs: critical temperatures [K]
                - Pcs: critical pressures [bar]

        Example
        -------
        .. code-block:: python

            import numpy as np

            import matplotlib.pyplot as plt

            from yaeos import PengRobinson76


            tc = np.array([369.83, 507.6])       # critical temperatures [K]
            pc = np.array([42.48, 30.25])        # critical pressures [bar]
            w = np.array([0.152291, 0.301261])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Two phase envelope calculation and plot
            env = model.phase_envelope_pt(
                np.array([0.5, 0.5]),
                t0=150.0,
                p0=1.0
            )

            plt.plot(env["Ts"], env["Ps"])
            plt.scatter(env["Tcs"], env["Pcs"])
        """
        ts, ps, tcs, pcs, xs, ys, kinds = yaeos_c.pt2_phase_envelope(
            self.id, z, kind=kind, t0=t0, p0=p0, max_points=max_points
        )

        msk = ~np.isnan(ts)
        msk_cp = ~np.isnan(tcs)

        res = {
            "T": ts[msk],
            "P": ps[msk],
            "Tc": tcs[msk_cp],
            "Pc": pcs[msk_cp],
            "x": xs[msk],
            "y": ys[msk],
            "kinds": kinds[msk],
        }

        return res

    def phase_envelope_px(
        self,
        z0,
        zi,
        temperature,
        kind="bubble",
        max_points=300,
        p0=10.0,
        a0=0.001,
        ns0=None,
        ds0=0.1,
    ):
        """Two phase envelope calculation (PX).
        Calculation of a phase envelope that starts at a given composition and
        its related to another composition with some proportion.

        Parameters
        ----------
        z0 : array_like
            Initial global mole fractions
        zi : array_like
            Final global mole fractions
        temperature : float
            Temperature [K]
        kind : str, optional
            Kind of saturation point to start the envelope calculation,
            defaults to "bubble". Options are
            - "bubble"
            - "dew"
        max_points : int, optional
            Envelope's maximum points to calculate (P, X), by default 300
        p0 : float, optional
            Initial guess for pressure [bar] for the saturation point of kind:
            `kind`, by default 10.0
        a0 : float, optional
            Initial molar fraction of composition `zi`, by default 0.001
        ns0 : int, optional
            Initial specified variable number, by default None.
            The the first `n=len(z)` values correspond to the K-values, where
            the last two values are the pressure and alpha.
        ds0 : float, optional
            Step for a, by default 0.1
        """

        if ns0 is None:
            ns0 = len(z0) + 2

        a, ps, xs, ys, acs, pcs, kinds = yaeos_c.px2_phase_envelope(
            self.id,
            z0=z0,
            zi=zi,
            kind=kind,
            max_points=max_points,
            p0=p0,
            a0=a0,
            t=temperature,
            ns0=ns0,
            ds0=ds0,
        )

        msk = ~np.isnan(ps)
        msk_cp = ~np.isnan(pcs)

        return {
            "a": a[msk],
            "P": ps[msk],
            "x": xs[msk],
            "y": ys[msk],
            "ac": acs[msk_cp],
            "Pc": pcs[msk_cp],
            "kind": kinds[msk],
        }

    def phase_envelope_tx(
        self,
        z0,
        zi,
        pressure,
        kind="bubble",
        max_points=300,
        t0=150.0,
        a0=0.001,
        ns0=None,
        ds0=0.1,
    ):
        """Two phase envelope calculation (TX).
        Calculation of a phase envelope that starts at a given composition and
        its related to another composition with some proportion.

        Parameters
        ----------
        z0 : array_like
            Initial global mole fractions
        zi : array_like
            Final global mole fractions
        pressure : float
            Pressure [K]
        kind : str, optional
            Kind of saturation point to start the envelope calculation,
            defaults to "bubble". Options are
            - "bubble"
            - "dew"
        max_points : int, optional
            Envelope's maximum points to calculate (P, X), by default 300
        t0 : float, optional
            Initial guess for temperature [K] for the saturation point of kind:
            `kind`, by default 150.0
        a0 : float, optional
            Initial molar fraction of composition `zi`, by default 0.001
        ns0 : int, optional
            Initial specified variable number, by default None.
            The the first `n=len(z)` values correspond to the K-values, where
            the last two values are the temperature and alpha.
        ds0 : float, optional
            Step for a, by default 0.1
        """

        if ns0 is None:
            ns0 = len(z0) + 2
        a, ts, xs, ys, acs, pcs, kinds = yaeos_c.tx2_phase_envelope(
            self.id,
            z0=z0,
            zi=zi,
            kind=kind,
            max_points=max_points,
            t0=t0,
            a0=a0,
            p=pressure,
            ns0=ns0,
            ds0=ds0,
        )

        msk = ~np.isnan(ts)
        msk_cp = ~np.isnan(pcs)

        return {
            "a": a[msk],
            "T": ts[msk],
            "x": xs[msk],
            "y": ys[msk],
            "ac": acs[msk_cp],
            "Pc": pcs[msk_cp],
            "kind": kinds[msk],
        }

    def phase_envelope_pt3(
        self,
        z,
        x0,
        y0,
        w0,
        beta0,
        t0,
        p0,
        specified_variable=None,
        first_step=None,
        max_points=1000,
    ):
        """
        Three-phase envelope tracing method.
        Calculation of a three-phase envelope that starts with an estimated
        compositions, pressure, temperature and phase fractions.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        x0 : array_like
            Initial phase x mole fractions
        y0 : array_like
            Initial phase y mole fractions
        w0 : array_like
            Initial incipient phase w mole fractions
        beta0 : float
            Initial phase fraction between x and y
        t0 : float
            Initial temperature [K]
        p0 : float
            Initial pressure [bar]
        specified_variable : int, optional
            Initial specified variable number, by default 2*len(z)+2
            (temperature).  The the first `n=(1,len(z))` values correspond to
            the K-values between phase x and w, the next `n=(len(z)+1,
            2*len(z))` are the K-values between phase y and w.  The last three
            values are pressure, temperature and beta.
        first_step : float, optional
            Step for the specified variable, by default 0.1
        max_points : int, optional
            Maximum number of points to calculate, by default 1000
        """

        if specified_variable is None:
            specified_variable = 2 * len(z) + 2

        if first_step is None:
            first_step = 0.1

        x, y, w, p, t, beta = yaeos_c.pt3_phase_envelope(
            self.id,
            z=z,
            x0=x0,
            y0=y0,
            w0=w0,
            beta0=beta0,
            t0=t0,
            p0=p0,
            max_points=max_points,
            ns0=specified_variable,
            ds0=first_step,
        )

        msk = ~np.isnan(t)

        return {
            "x": x[msk],
            "y": y[msk],
            "w": w[msk],
            "P": p[msk],
            "T": t[msk],
            "beta": beta[msk],
        }

    # ==============================================================
    # Stability analysis
    # --------------------------------------------------------------
    def stability_analysis(self, z, pressure, temperature):
        """Perform stability analysis.
        Find all the possible minima values that the :math:`tm` function,
        defined by Michelsen and Mollerup.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        pressure : float
            Pressure [bar]
        temperature : float
            Temperature [K]

        Returns
        -------
        dict
            Stability analysis result dictionary with keys:
            - w: value of the test phase that minimizes the :math:`tm` function.
            - tm: minimum value of the :math:`tm` function.
        dict
            All found minimum values of the :math:`tm` function and the
            corresponding test phase mole fractions.
            - w: all values of :math:`w` that minimize the :math:`tm` function
            - tm: all values found minima of the :math:`tm` function"""
        (w_min, tm_min, all_mins, all_mins_w) = yaeos_c.stability_zpt(
            id=self.id, z=z, p=pressure, t=temperature
        )

        return {"w": w_min, "tm": tm_min}, {
            "tm": all_mins,
            "w": all_mins_w,
        }

    def stability_tm(self, z, w, pressure, temperature):
        """Calculate the :math:`tm` function.

        Calculate the :math:`tm` function, defined by Michelsen and Mollerup.
        If this value is negative, it means that the feed with composition `z`
        is unstable.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        w : array_like
            Test Phase mole fractions
        pressure : float
            Pressure [bar]
        temperature : float
            Temperature [K]

        Returns
        -------
        float
            Value of the :math:`tm` function
        """
        return yaeos_c.tm(id=self.id, z=z, w=w, p=pressure, t=temperature)

    # ==============================================================
    # Critical points and lines
    # --------------------------------------------------------------
    def critical_point(self, z0, zi=[0, 0], ns=1, s=0, max_iters=100) -> dict:
        """Critical point calculation.

        Calculate the critical point of a mixture. At a given composition.

        Parameters
        ----------
        z0: array_like
            Mole fractions of original fluid
        zi: array_like
            Mole fractinos of other fluid
        ns: int
            Number of specification
        S: float
            Specification value
        max_iters: int, optional

        Returns
        -------
        dict
            Critical point calculation result dictionary with keys:
                - Tc: critical temperature [K]
                - Pc: critical pressure [bar]
                - Vc: critical volume [L]
        """
        *x, t, p, v = yaeos_c.critical_point(
            self.id, z0=z0, zi=zi, spec=ns, s=s, max_iters=max_iters
        )

        return {"x": x, "Tc": t, "Pc": p, "Vc": v}

    def critical_line(
        self,
        z0,
        zi,
        ns=1,
        s=1e-5,
        ds0=1e-2,
        a0=1e-5,
        max_points=1000,
        stop_pressure=2500,
    ):
        """Critical Line calculation.
        Calculate the critical line between two compositions

        Parameters
        ----------
        z0: array_like
            Initial global mole fractions
        zi: array_like
            Final global mole fractions
        ns: int, optional
            Specified variable number, by default 1
        s: float, optional
            Specified value, by default 1e-5
        ds0: float, optional
            Step for molar fraction of composition `i`
        a0: float, optional
            Initial molar fraction of composition `i`
        max_points: int, optional
            Maximum number of points to calculate
        stop_pressure: float, optional
            Stop when reaching this pressure value
        """

        alphas, vs, ts, ps = yaeos_c.critical_line(
            self.id,
            ns=ns,
            ds0=ds0,
            a0=a0,
            s=s,
            z0=z0,
            zi=zi,
            max_points=max_points,
            stop_pressure=stop_pressure,
        )

        msk = ~np.isnan(ts)

        return {"a": alphas[msk], "T": ts[msk], "P": ps[msk], "V": vs[msk]}

    def __del__(self) -> None:
        """Delete the model from the available models list (Fortran side)."""
        yaeos_c.make_available_ar_models_list(self.id)
