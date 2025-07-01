"""yaeos Python API core module.

ArModel and GeModel abstract classes definition. Also, the implementation of
the models' thermoprops methods.
"""

from abc import ABC
from typing import Union

from intersect import intersection

import numpy as np

from yaeos.lib import yaeos_c

from yaeos.envelopes import PTEnvelope, PXEnvelope, TXEnvelope

from yaeos.constants import root_kinds

from warnings import warn


MAX_POINTS_ENVELOPES = 1000000


def adjust_root_kind(number_of_phases, kinds_x=None, kind_w=None):
    """Convert the the kinds of each phase to the corresponding value.

    The C interface of `yaeos` expects the kinds of each phase to be defined
    as integer values, so this function converts the kinds of each phase
    to the corresponding integer value. If the kind is not specified, it
    defaults to "stable" for all phases.

    Parameters
    ----------
    number_of_phases : int
        Number of phases in the system, besides de reference phase
    kinds_x : list, optional
        Kinds of the phases in the system, by default None
    kind_w : str, optional
        Kind of the test phase, by default None
    Returns
    -------
    tuple
        kinds_x_out : list
            List of kinds for each phase in the system
        kind_w_out : str
            Kind of the test phase
    """
    if kinds_x:
        kinds_x_out = [root_kinds[kind] for kind in kinds_x]
    else:
        kinds_x_out = [root_kinds["stable"] for _ in range(number_of_phases)]

    if kind_w:
        kind_w_out = root_kinds[kind_w]
    else:
        kind_w_out = root_kinds["stable"]

    return kinds_x_out, kind_w_out


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

        possible_derivatives = [dt, dt2, dn, dtn]
        all_none = all([d is None for d in possible_derivatives])

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

    def stability_analysis(self, z, temperature):
        """Perform stability analysis.

        Find all the possible minima values that the :math:`tm` function,
        defined by Michelsen and Mollerup.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        temperature : float
            Temperature [K]

        Returns
        -------
        dict
            Stability analysis result dictionary with keys:
            - w: value of the test phase that minimizes the :math:`tm` function
            - tm: minimum value of the :math:`tm` function.
        dict
            All found minimum values of the :math:`tm` function and the
            corresponding test phase mole fractions.
            - w: all values of :math:`w` that minimize the :math:`tm` function
            - tm: all values found minima of the :math:`tm` function"""
        (w_min, tm_min, all_mins) = yaeos_c.stability_zt_ge(
            id=self.id, z=z, t=temperature
        )

        all_mins_w = all_mins[:, : len(z)]
        all_mins = all_mins[:, -1]

        return {"w": w_min, "tm": tm_min}, {"tm": all_mins, "w": all_mins_w}

    def flash_t(self, z, temperature: float, k0=None) -> dict:
        """Two-phase split with specification of temperature and pressure.

        Parameters
        ----------
        z : array_like
            Global mole fractions
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
                - T: temperature [K]
                - beta: light phase fraction
        """
        if k0 is None:
            mintpd, _ = self.stability_analysis(z, temperature)
            k0 = mintpd["w"] / np.array(z)

        x, y, pressure, temperature, volume_x, volume_y, beta = (
            yaeos_c.flash_ge(self.id, z, t=temperature, k0=k0)
        )

        flash_result = {
            "x": x,
            "y": y,
            "Vx": volume_x,
            "Vy": volume_y,
            "T": temperature,
            "beta": beta,
        }

        return flash_result

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
            self.id, component, stop_p=stop_pressure, stop_t=stop_temperature
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
                - beta: light phase fraction. If beta is -1 flash was not
                successful.

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

    def flash_vt(self, z, volume: float, temperature: float, k0=None) -> dict:
        """Two-phase split with specification of temperature and volume.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        volume : float
            Molar volume [L/mol]
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
                - Vx: heavy phase molar volume [L/mol]
                - Vy: light phase molar volume [L/mol]
                - P: pressure [bar]
                - T: temperature [K]
                - beta: light phase fraction. If beta is -1 flash was not
                successful.

        Example
        -------
        .. code-block:: python

            import numpy as np

            from yaeos import PengRobinson76


            tc = np.array([507.6, 658.0])        # critical temperatures [K]
            pc = np.array([30.25, 18.20])        # critical pressures [bar]
            w = np.array([0.301261, 0.576385])   # acentric factors

            model = PengRobinson76(tc, pc, w)

            # Flash calculation
            # will print:
            # {
            #     'x': array([0.26308567, 0.73691433]),
            #     'y': array([0.95858707, 0.04141293]),
            #     'Vx': 0.2417828483590114,
            #     'Vy': 31.706890870110417,
            #     'P': 1.0001131874567775,
            #     'T': 393.15,
            #     'beta': 0.34063818069513246
            # }

            print(model.flash_vt([0.5, 0.5], 10.96, 393.15))
        """
        if k0 is None:
            k0 = [0 for i in range(len(z))]

        x, y, pressure, temperature, volume_x, volume_y, beta = (
            yaeos_c.flash_vt(self.id, z, v=volume, t=temperature, k0=k0)
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
        self,
        z,
        temperature: float,
        kind: str = "bubble",
        p0: float = 0,
        y0=None,
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
        p0: float, optional
            Initial guess for pressure [bar]
        y0: array_like, optional
            Initial guess for the incipient phase, by default None
            (will use k_wilson correlation)

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
        if y0 is None:
            y0 = np.zeros_like(z)

        p, x, y, volume_x, volume_y, beta = yaeos_c.saturation_pressure(
            id=self.id, z=z, t=temperature, kind=kind, p0=p0, y0=y0
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
        self, z, pressure: float, kind: str = "bubble", t0: float = 0, y0=None
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
        t0: float, optional
            Initial guess for temperature [K]
        y0: array_like, optional
            Initial guess for the incipient phase, by default None
            (will use k_wilson correlation)

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
        if y0 is None:
            y0 = np.zeros_like(z)

        t, x, y, volume_x, volume_y, beta = yaeos_c.saturation_temperature(
            id=self.id, z=z, p=pressure, kind=kind, t0=t0, y0=y0
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

    # =========================================================================
    # Phase envelopes
    # -------------------------------------------------------------------------
    def phase_envelope_pt(
        self,
        z,
        kind: str = "bubble",
        max_points: int = MAX_POINTS_ENVELOPES,
        t0: float = 150.0,
        p0: float = 1.0,
        w0=None,
        stop_pressure: float = 2500,
        ds0: float = 0.001,
    ) -> PTEnvelope:
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
            - "liquid-liquid"
        max_points : int, optional
            Envelope's maximum points to calculate (T, P)
        t0 : float, optional
            Initial guess for temperature [K] for the saturation point of kind:
            `kind`, by default 150.0
        p0 : float, optional
            Initial guess for pressure [bar] for the saturation point of kind:
            `kind`, by default 1.0
        w0 : array_like, optional
            Initial guess for the incipient phase mole fractions,
            by default None In the case of bubble and dew line calculations, it
            will use the k_wilson correlation. In the case of liquid-liquid
            envelope it will make a search for the first unstable component
            when decreasing temperature at the given pressure.
        stop_pressure : float, optional
            Stop on pressures above stop_pressure [bar], by default 2500.0.
            If the the initial guess pressure is above this value, the
            calculation will stop immediately.
        ds0: float, optional
            Step for the first specified variable, by default 0.001. The
            specified variable is the temperature for bubble and dew lines, and
            pressure for liquid-liquid lines. For bubble and dew lines, the
            step is positive, while for liquid-liquid lines it is negative.

        Returns
        -------
        PTEnvelope
            PTEnvelope object with the phase envelope information.

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

            plt.plot(env["T"], env["P"])
            plt.scatter(env["Tc"], env["Pc"])
        """

        ds0 = 0.001

        if kind == "bubble":
            sat = self.saturation_pressure(z, t0, kind=kind, p0=p0)
            w0 = sat["y"]
            ns0 = len(z) + 3
            p0 = sat["P"]
            kinds_x = ["liquid"]
            kind_w = "vapor"
        elif kind == "dew":
            sat = self.saturation_temperature(z, p0, kind=kind, t0=t0)
            w0 = sat["x"]
            ns0 = len(z) + 3
            t0 = sat["T"]
            kinds_x = ["vapor"]
            kind_w = "liquid"
        elif kind == "liquid-liquid":
            ns0 = len(z) + 2
            if w0 is None:
                # =============================================================
                # Find an initialization for the liquid-liquid envelope
                # -------------------------------------------------------------
                ts = []
                for i in range(len(z)):
                    w0 = np.zeros_like(z)
                    w0 += 1e-5
                    w0[i] = 1 - np.sum(w0[1:])
                    for t in np.linspace(1000, 100, 25):
                        tm = self.stability_tm(z, w0, p0, t)
                        if tm < -0.01:
                            ts.append(t)
                            break
                if len(ts) == 0:
                    warn("No liquid-liquid region found.")
                    return None
                i = np.argmin(ts)
                t0 = ts[i]
                w0 = np.zeros_like(z)
                w0 += 1e-5
                w0[i] = 1 - np.sum(w0[1:])

            kinds_x = ["liquid"]
            kind_w = "liquid"
            ds0 = -ds0

        envelope = self.phase_envelope_pt_mp(
            z,
            x_l0=[z],
            w0=w0,
            betas0=[1],
            t0=t0,
            p0=p0,
            ns0=ns0,
            ds0=ds0,
            max_points=max_points,
            stop_pressure=stop_pressure,
            kinds_x=kinds_x,
            kind_w=kind_w,
        )

        return envelope

    def phase_envelope_px(
        self,
        z0,
        zi,
        temperature,
        kind="bubble",
        max_points=MAX_POINTS_ENVELOPES,
        p0=10.0,
        w0=None,
        a0=1e-2,
        ns0=None,
        ds0=1e-5,
    ) -> PXEnvelope:
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
            Envelope's maximum points to calculate.
        p0 : float, optional
            Initial guess for pressure [bar] for the saturation point of kind:
            `kind`, by default 10.0
        a0 : float, optional
            Initial molar fraction of composition `zi`, by default 0.001
        ns0 : int, optional
            Initial specified variable number, by default None.
            The the first `n=len(z)` values correspond to the K-values,
            `len(z)+1` is the main phase molar fraction (ussually one) and
            the last two values are the pressure and alpha.
        ds0 : float, optional
            Step for the first specified variable, by default 0.01
        """
        if ns0 is None:
            ns0 = len(z0) + 3

        zi = np.array(zi)
        z0 = np.array(z0)

        if w0 is None:
            w0 = np.zeros_like(z0)

        z = a0 * zi + (1 - a0) * z0
        sat = self.saturation_pressure(
            z, temperature=temperature, kind=kind, p0=p0, y0=w0
        )

        if kind == "dew":
            w0 = sat["x"]
            kind_x = ["vapor"]
            kind_w = "liquid"
        else:
            w0 = sat["y"]
            kind_x = ["liquid"]
            kind_w = "vapor"

        p0 = sat["P"]

        envelope = self.phase_envelope_px_mp(
            z0,
            zi,
            temperature,
            x_l0=[z],
            w0=w0,
            betas0=[1],
            p0=p0,
            alpha0=a0,
            ns0=ns0,
            ds0=ds0,
            max_points=max_points,
            kinds_x=kind_x,
            kind_w=kind_w,
        )

        return envelope

    def phase_envelope_tx(
        self,
        z0,
        zi,
        pressure,
        kind="bubble",
        max_points=MAX_POINTS_ENVELOPES,
        t0=150.0,
        a0=0.001,
        ns0=None,
        ds0=0.1,
        w0=None,
    ) -> TXEnvelope:
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
            Pressure [bar]
        kind : str, optional
            Kind of saturation point to start the envelope calculation,
            defaults to "bubble". Options are
            - "bubble"
            - "dew"
        max_points : int, optional
            Envelope's maximum points to calculate (P, X).
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
        w0 : array_like, optional
            Initial guess for the incipient phase, by default None
            (will use k_wilson correlation)
        """
        zi = np.array(zi)
        z0 = np.array(z0)

        if not ns0:
            ns0 = len(z0) + 3

        if w0 is None:
            w0 = np.zeros_like(z0)

        z = a0 * zi + (1 - a0) * z0
        sat = self.saturation_temperature(
            z, pressure=pressure, kind=kind, t0=t0, y0=w0
        )

        if kind == "dew":
            w0 = sat["x"]
            kind_x = ["vapor"]
            kind_w = "liquid"
        else:
            w0 = sat["y"]
            kind_x = ["liquid"]
            kind_w = "vapor"

        envelope = self.phase_envelope_tx_mp(
            z0=z0,
            zi=zi,
            p=pressure,
            x_l0=[z],
            w0=w0,
            betas0=[1],
            t0=t0,
            alpha0=a0,
            ns0=ns0,
            ds0=ds0,
            max_points=max_points,
            kinds_x=kind_x,
            kind_w=kind_w,
        )

        return envelope

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
        kinds_x=None,
        kind_w=None,
        max_points=MAX_POINTS_ENVELOPES,
        stop_pressure=2500,
    ) -> PTEnvelope:
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
        kinds_x : list, optional
            Kinds of the main phases, by default None (will use "stable")
        kind_w : str, optional
            Kind of the reference phase, by default None (will use "stable")
        max_points : int, optional
            Maximum number of points to calculate
        stop_pressure : float, optional
            Stop at pressure above stop_pressure [bar], default 2500
        """

        np = 2
        if specified_variable is None:
            specified_variable = 2 * len(z) + np + 2

        if first_step is None:
            first_step = 0.1

        envelope = self.phase_envelope_pt_mp(
            z=z,
            x_l0=[x0, y0],
            w0=w0,
            betas0=[1 - beta0, beta0],
            t0=t0,
            p0=p0,
            ns0=specified_variable,
            ds0=first_step,
            beta_w=0,
            kinds_x=kinds_x,
            kind_w=kind_w,
            max_points=max_points,
            stop_pressure=stop_pressure,
        )

        return envelope

    def phase_envelope_px3(
        self,
        z0,
        zi,
        T,
        x0,
        y0,
        w0,
        beta0,
        a0,
        p0,
        specified_variable=None,
        first_step=None,
        max_points=MAX_POINTS_ENVELOPES,
        kinds_x=None,
        kind_w=None,
    ) -> PXEnvelope:
        """
        Three-phase envelope tracing method.

        Calculation of a three-phase envelope that starts with an estimated
        compositions, pressure, temperature and phase fractions.

        Parameters
        ----------
        z0 : array_like
            Global mole fractions of the original fluid
        zi : array_like
            Global mole fractions of the other fluid
        x0 : array_like
            Initial phase x mole fractions
        y0 : array_like
            Initial phase y mole fractions
        w0 : array_like
            Initial incipient phase w mole fractions
        beta0 : float
            Initial phase fraction between x and y
        a0 : float
            Initial molar fraction of the other fluid
        p0 : float
            Initial pressure [bar]
        specified_variable : int, optional
            Initial specified variable number, by default 2*len(z)+2
            (temperature).  The the first `n=(1,len(z))` values correspond to
            the K-values between phase x and w, the next `n=(len(z)+1,
            2*len(z))` are the K-values between phase y and w.  The last three
            values are pressure, a and beta.
        first_step : float, optional
            Step for the specified variable, by default 0.1
        max_points : int, optional
            Maximum number of points to calculate
        kinds_x : list, optional
            Kinds of the main phases, by default None (will use "stable")
            options can be - "stable", "liquid", "vapor"
        kind_w : str, optional
            Kind of the reference phase, by default None (will use "stable")
            options can be - "stable", "liquid", "vapor"
        """
        if specified_variable is None:
            specified_variable = 2 * len(z0) + 2

        if first_step is None:
            first_step = 0.1

        kinds_x, kind_w = adjust_root_kind(
            number_of_phases=2, kinds_x=kinds_x, kind_w=kind_w
        )

        envelope = self.phase_envelope_px_mp(
            z0=z0,
            zi=zi,
            t=T,
            x_l0=[x0, y0],
            w0=w0,
            betas0=[1 - beta0, beta0],
            p0=p0,
            alpha0=a0,
            ns0=specified_variable,
            ds0=first_step,
            max_points=max_points,
            kinds_x=kinds_x,
            kind_w=kind_w,
        )
        return envelope

    def phase_envelope_pt_mp(
        self,
        z,
        x_l0,
        w0,
        betas0,
        p0,
        t0,
        ns0,
        ds0,
        beta_w=0,
        kinds_x=None,
        kind_w=None,
        max_points=MAX_POINTS_ENVELOPES,
        stop_pressure=1000,
    ) -> PTEnvelope:
        """Multi-phase envelope."""
        x_l0 = np.array(x_l0)

        number_of_phases = x_l0.shape[0]

        kinds_x, kind_w = adjust_root_kind(
            number_of_phases=number_of_phases, kinds_x=kinds_x, kind_w=kind_w
        )

        x_ls, ws, betas, ps, ts, iters, ns, x_kinds, w_kinds, pcs, tcs = (
            yaeos_c.pt_mp_phase_envelope(
                id=self.id,
                np=number_of_phases,
                z=z,
                x_l0=x_l0,
                w0=w0,
                betas0=betas0,
                t0=t0,
                p0=p0,
                ns0=ns0,
                ds0=ds0,
                beta_w=beta_w,
                kinds_x=kinds_x,
                kind_w=kind_w,
                max_points=max_points,
                stop_pressure=stop_pressure,
            )
        )

        x_kinds = x_kinds.astype(str)
        w_kinds = w_kinds.astype(str)
        x_kinds = np.char.strip(x_kinds)
        w_kinds = np.char.strip(w_kinds)

        envelope = PTEnvelope(
            global_composition=z,
            main_phases_compositions=x_ls,
            reference_phase_compositions=ws,
            reference_phase_kinds=w_kinds,
            main_phases_kinds=x_kinds,
            main_phases_molar_fractions=betas,
            pressures=ps,
            temperatures=ts,
            iterations=iters,
            specified_variable=ns,
            critical_pressures=pcs,
            critical_temperatures=tcs,
        )

        return envelope

    def phase_envelope_px_mp(
        self,
        z0,
        zi,
        t,
        x_l0,
        w0,
        betas0,
        p0,
        ns0,
        ds0,
        alpha0=0,
        beta_w=0,
        max_points=MAX_POINTS_ENVELOPES,
        kinds_x=None,
        kind_w=None,
    ) -> PXEnvelope:
        """Multi-phase PX envelope.

        Calculate a phase envelope with a preselected ammount of phases.

        Parameters
        ----------
        z0: float, array_like
            Original Fluid.
        zi: float, array_like
            Other fluid.
        t: float
            Temperature [K]
        x_l0: float, matrix [number of phases, number of components]
            A matrix where each row is the composition of a main phase.
            Guess for first point.
        w0: flot, array_like
            Composition of the reference (ussually incipient) phase.
            Guess for first point
        betas0: float, array_like
            Molar fraction of each main phase. Guess for first point
        p0: float
            Pressure guess for first point [bar]
        ns0: int
            Initial variable to specifiy.
            From 1 to `(number_of_phases*number_of_components)` corresponds to
            each composition.
            From `number_of_phases*number_of_components` to
            `number_of_phases*number_of_components + number_of_phases`
            corresponds to each beta value of the main phases.
            The last two posibilities are the pressure and molar relation
            between the two fluids, respectively.
        kinds_x: list(str), optional
            List of kinds of main phases, defaults to stable. options are:
            - "stable"
            - "liquid"
            - "vapor"
        kinds_w: list(str), optional
            Kind of reference phase, defaults to stable. options are:
            - "stable"
            - "liquid"
            - "vapor"
        """
        x_l0 = np.array(x_l0, order="F")

        number_of_phases = x_l0.shape[0]

        kinds_x, kind_w = adjust_root_kind(
            number_of_phases=number_of_phases, kinds_x=kinds_x, kind_w=kind_w
        )

        x_ls, ws, betas, ps, alphas, iters, ns, x_kinds, w_kinds, pcs, acs = (
            yaeos_c.px_mp_phase_envelope(
                id=self.id,
                np=number_of_phases,
                z0=z0,
                zi=zi,
                x_l0=x_l0,
                w0=w0,
                betas0=betas0,
                t=t,
                beta_w=beta_w,
                kinds_x=kinds_x,
                kind_w=kind_w,
                alpha0=alpha0,
                p0=p0,
                ns0=ns0,
                ds0=ds0,
                max_points=max_points,
            )
        )

        return PXEnvelope(
            temperature=t,
            global_composition_0=z0,
            global_composition_i=zi,
            main_phases_compositions=x_ls,
            reference_phase_compositions=ws,
            main_phases_molar_fractions=betas,
            pressures=ps,
            alphas=alphas,
            iterations=iters,
            specified_variable=ns,
            critical_pressures=pcs,
            critical_alphas=acs,
            main_phases_kinds=x_kinds,
            reference_phase_kinds=w_kinds,
        )

    def phase_envelope_tx_mp(
        self,
        z0,
        zi,
        p,
        x_l0,
        w0,
        betas0,
        t0,
        ns0,
        ds0,
        alpha0=0,
        beta_w=0,
        max_points=MAX_POINTS_ENVELOPES,
        kinds_x=None,
        kind_w=None,
    ) -> TXEnvelope:
        """Multi-phase envelope."""

        x_l0 = np.array(x_l0, order="F")

        number_of_phases = x_l0.shape[0]

        kinds_x, kind_w = adjust_root_kind(
            number_of_phases=number_of_phases, kinds_x=kinds_x, kind_w=kind_w
        )

        (
            x_ls,
            ws,
            betas,
            ts,
            alphas,
            iters,
            ns,
            main_kinds,
            ref_kinds,
            tcs,
            acs,
        ) = yaeos_c.tx_mp_phase_envelope(
            id=self.id,
            np=number_of_phases,
            z0=z0,
            zi=zi,
            p=p,
            beta_w=beta_w,
            kinds_x=kinds_x,
            kind_w=kind_w,
            x_l0=x_l0,
            w0=w0,
            betas0=betas0,
            alpha0=alpha0,
            t0=t0,
            ns0=ns0,
            ds0=ds0,
            max_points=max_points,
        )

        return TXEnvelope(
            pressure=p,
            global_composition_0=z0,
            global_composition_i=zi,
            main_phases_compositions=x_ls,
            reference_phase_compositions=ws,
            main_phases_molar_fractions=betas,
            temperatures=ts,
            alphas=alphas,
            iterations=iters,
            specified_variable=ns,
            critical_temperatures=tcs,
            critical_alphas=acs,
            main_phases_kinds=main_kinds,
            reference_phase_kinds=ref_kinds,
        )

    def phase_envelope_pt_from_dsp(
        self,
        z,
        env1: PTEnvelope,
        env2: PTEnvelope,
        dbeta0=1e-5,
        max_points=MAX_POINTS_ENVELOPES,
    ) -> list:
        """Calculate PT phase envelopes from a DSP.

        This method calculates the phase envelope at the intersection of two
        PT envelopes, `env1` and `env2`.

        Parameters
        ----------
        z : array_like
            Global mole fractions
        env1 : PTEnvelope
            First PT envelope object
        env2 : PTEnvelope
            Second PT envelope object
        dbeta0 : float, optional
            initial step for the beta values, by default 1e-5
        max_points : int, optional
            Maximum number of points to calculate
        """
        nc = env1.number_of_components
        phases = env1.number_of_phases + 1

        Ts, Ps = intersection(env1["T"], env1["P"], env2["T"], env2["P"])

        dsps = []
        for Tdsp, Pdsp in zip(Ts, Ps):
            env1_loc = np.argmin(
                np.abs(env1["T"] - Tdsp) + np.abs(env1["P"] - Pdsp)
            )
            env2_loc = np.argmin(
                np.abs(env2["T"] - Tdsp) + np.abs(env2["P"] - Pdsp)
            )

            betas_1 = env1.main_phases_molar_fractions[env1_loc, :]
            betas_2 = env2.main_phases_molar_fractions[env2_loc, :]

            w0 = env1.reference_phase_compositions[env1_loc]
            y0 = env2.reference_phase_compositions[env2_loc]

            x_l1 = np.vstack(
                (env1.main_phases_compositions[env1_loc, :, :], y0)
            )
            x_l2 = np.vstack(
                (env1.main_phases_compositions[env1_loc, :, :], w0)
            )

            # Convert the kinds to the correct format. Saving first as *_ to
            # avoid issues.
            kinds_x_1 = env1.main_phases_kinds[env1_loc, :]
            kinds_x_2 = env2.main_phases_kinds[env2_loc, :]
            kind_w_1 = env1.reference_phase_kinds[env1_loc]
            kind_w_2 = env2.reference_phase_kinds[env2_loc]
            dsp_1 = self.phase_envelope_pt_mp(
                z=z,
                x_l0=x_l1,
                w0=w0,
                betas0=[*betas_1, 0],
                p0=Pdsp,
                t0=Tdsp,
                ns0=phases * nc + phases,
                ds0=dbeta0,
                beta_w=0,
                kinds_x=[*kinds_x_1, kind_w_1],
                kind_w=kind_w_2,
                max_points=max_points,
            )

            dsp_2 = self.phase_envelope_pt_mp(
                z=z,
                x_l0=x_l2,
                w0=y0,
                betas0=[*betas_2, 0],
                p0=Pdsp,
                t0=Tdsp,
                ns0=phases * nc + phases,
                ds0=dbeta0,
                beta_w=0,
                kinds_x=[*kinds_x_2, kind_w_2],
                kind_w=kind_w_1,
                max_points=max_points,
            )

            dsps.append([dsp_1, dsp_2])

        return dsps

    def phase_envelope_px_from_dsp(
        self, z0, zi, env1: PXEnvelope, env2: PXEnvelope, dbeta0=1e-5
    ) -> list:
        """Calculate PX phase envelopes from a DSP.

        This method calculates the phase envelope at the intersection of two
        PX envelopes, `env1` and `env2`.
        Parameters
        ----------
        z0 : array_like
            Global mole fractions of the original fluid
        zi : array_like
            Global mole fractions of the other fluid
        env1 : PXEnvelope
            First PX envelope object
        env2 : PXEnvelope
            Second PX envelope object
        dbeta0 : float, optional
            Initial step for the beta values, by default 1e-5
        Returns
        -------
        list
            List of lists of two PXEnvelope objects, one for each intersection
            point.
        """

        nc = env1.number_of_components
        phases = env1.number_of_phases + 1
        dsps = intersection(env1["a"], env1["P"], env2["a"], env2["P"])
        dsps = []

        for adsp, Pdsp in zip(dsps[0], dsps[1]):
            env1_loc = np.argmin(abs(env1["a"] - adsp) + abs(env1["P"]) - Pdsp)
            env2_loc = np.argmin(abs(env2["a"] - adsp) + abs(env2["P"]) - Pdsp)

            betas_1 = env1.main_phases_molar_fractions[env1_loc, :]
            betas_2 = env2.main_phases_molar_fractions[env2_loc, :]

            w0 = env1.reference_phase_compositions[env1_loc]
            y0 = env2.reference_phase_compositions[env2_loc]

            x_l1 = np.vstack(
                (env1.main_phases_compositions[env1_loc, :, :], y0)
            )
            x_l2 = np.vstack(
                (env1.main_phases_compositions[env1_loc, :, :], w0)
            )

            dsp_1 = self.phase_envelope_px_mp(
                z0=z0,
                zi=zi,
                t=T,
                x_l0=x_l1,
                w0=w0,
                betas0=[*betas_1, 0],
                p0=Pdsp,
                alpha0=adsp,
                ns0=phases * nc + phases,
                ds0=dbeta0,
                beta_w=0,
                max_points=MAX_POINTS_ENVELOPES,
            )

            dsp_2 = self.phase_envelope_px_mp(
                z0=z0,
                zi=zi,
                t=T,
                x_l0=x_l2,
                w0=y0,
                betas0=[*betas_2, 0],
                p0=Pdsp,
                alpha0=adsp,
                ns0=phases * nc + phases,
                ds0=dbeta0,
                beta_w=0,
                max_points=MAX_POINTS_ENVELOPES,
            )

            dsps.append([dsp_1, dsp_2])

        return dsps

    def isopleth(
        self,
        z,
        three_phase=True,
        dew_start=(500, 0.01),
        bubble_start=(200, 10),
        max_points=MAX_POINTS_ENVELOPES,
        delta_dew_2ph=0.01,
        delta_bub_2ph=0.01,
        delta_dsp_3ph=0.01,
        stop_pressure=2500,
    ):

        dew_point = self.saturation_temperature(
            z, pressure=dew_start[1], kind="dew", t0=dew_start[0]
        )
        bub_point = self.saturation_pressure(
            z, temperature=bubble_start[0], kind="bubble", p0=bubble_start[1]
        )

        dew_line = self.phase_envelope_pt_mp(
            z=z,
            x_l0=[z],
            w0=dew_point["x"],
            betas0=[1],
            p0=dew_point["P"],
            t0=dew_point["T"],
            ns0=len(z) + 2,
            ds0=delta_dew_2ph,
            max_points=max_points,
            stop_pressure=stop_pressure,
        )

        bub_line = self.phase_envelope_pt_mp(
            z=z,
            x_l0=[z],
            w0=bub_point["y"],
            betas0=[1],
            p0=bub_point["P"],
            t0=bub_point["T"],
            ns0=len(z) + 2,
            ds0=delta_bub_2ph,
            max_points=max_points,
            stop_pressure=stop_pressure,
        )

        liq = self.phase_envelope_pt(
            z, kind="liquid-liquid", t0=500, p0=2000, max_points=2
        )

        if len(liq["T"]) > 0:
            liq_line = self.phase_envelope_pt_mp(
                z=z,
                x_l0=[z],
                w0=liq.reference_phase_compositions[0],
                betas0=[1],
                p0=liq["P"][0],
                t0=liq["T"][0],
                ns0=len(z) + 2,
                ds0=-0.01,
                max_points=max_points,
                stop_pressure=1e10,
            )
        else:
            liq_line = None

        dsps_db = intersection(
            dew_line["T"],
            dew_line["P"],
            bub_line["T"],
            bub_line["P"],
        )

        if liq_line:
            dsps_dl = intersection(
                dew_line["T"],
                dew_line["P"],
                liq_line["T"],
                liq_line["P"],
            )

            dsps_bl = intersection(
                bub_line["T"],
                bub_line["P"],
                liq_line["T"],
                liq_line["P"],
            )

            dsps_set = {
                "dl": [dsps_dl, dew_line, liq_line],
                "db": [dsps_db, dew_line, bub_line],
                "bl": [dsps_bl, bub_line, liq_line],
            }
        else:
            dsps_set = {"db": [dsps_db, dew_line, bub_line]}

        dew_locs = []
        bub_locs = []
        liq_locs = []

        three_phase_envs = []
        stable_lines = {"3ph": [], "2ph": []}

        if three_phase:
            dew_line_stable = dew_line
            bub_line_stable = bub_line
            liq_line_stable = liq_line

            for dsp_name in dsps_set:
                # Order the DSPs wrt temperature
                dsps = dsps_set[dsp_name][0]
                env_1, env_2 = dsps_set[dsp_name][1], dsps_set[dsp_name][2]

                idx = dsps[0].argsort()

                dsps = (dsps[0][idx], dsps[1][idx])

                if 0 < len(dsps[0]) <= 2:
                    for temperature, pressure in zip(dsps[0], dsps[1]):
                        env_1_loc = np.argmin(
                            np.abs(env_1["T"] - temperature)
                            + np.abs(env_1["P"] - pressure)
                        )
                        env_2_loc = np.argmin(
                            np.abs(env_2["T"] - temperature)
                            + np.abs(env_2["P"] - pressure)
                        )

                        dew_locs.append(env_1_loc)
                        bub_locs.append(env_2_loc)

                        w_env_1 = env_1.reference_phase_compositions[env_1_loc]
                        w_env_2 = env_2.reference_phase_compositions[env_2_loc]

                        env1 = self.phase_envelope_pt_mp(
                            z,
                            x_l0=np.array([z, w_env_1]),
                            w0=w_env_2,
                            betas0=[1, 0],
                            t0=temperature,
                            p0=pressure,
                            ns0=2 * len(z) + 2,
                            ds0=delta_dsp_3ph,
                            max_points=max_points,
                            stop_pressure=max([pressure * 2, stop_pressure]),
                        )

                        env2 = self.phase_envelope_pt_mp(
                            z,
                            x_l0=np.array([z, w_env_2]),
                            w0=w_env_1,
                            betas0=[1, 0],
                            t0=temperature,
                            p0=pressure,
                            ns0=2 * len(z) + 2,
                            ds0=delta_dsp_3ph,
                            max_points=max_points,
                            stop_pressure=max([pressure * 2, stop_pressure]),
                        )

                        three_phase_envs.append((env1, env2))
                    # if len(dew_locs) == 2:
                    # msk = np.array([False] * len(dew_line))
                    # msk[: dew_locs[1] + 1] = True
                    # msk[dew_locs[0] :] = True
                    # dew_line_stable *= np.nan
                    # dew_line_stable[msk] = dew_line[msk]

                    # msk = np.array([False] * len(bub_line))
                    # msk[bub_locs[0] : bub_locs[1] + 1] = True
                    # bub_line_stable *= np.nan
                    # bub_line_stable[msk] = bub_line[msk]

                elif len(dsps[0]) == 0:

                    k0 = dew_line.reference_phase_compositions[0, :] / z
                    flash = self.flash_pt(
                        z,
                        pressure=bub_line["P"][0],
                        temperature=bub_line["T"][0],
                        k0=k0,
                    )

                    x_l0 = [z, flash["y"]]
                    w0 = bub_line.reference_phase_compositions[0, :]
                    betas = [1 - flash["beta"], flash["beta"]]

                    bubble_isolated = self.phase_envelope_pt_mp(
                        z=z,
                        x_l0=x_l0,
                        w0=w0,
                        betas0=betas,
                        p0=flash["P"],
                        t0=flash["T"],
                        ns0=2 * len(z) + 4,
                        ds0=delta_bub_2ph,
                        max_points=max_points,
                    )

                    x_l0 = [z, bub_line.reference_phase_compositions[0, :]]
                    w0 = flash["y"]
                    idx = np.argmax(w0)

                    flash = self.flash_pt(z, flash["P"], flash["T"])
                    betas = [1 - flash["beta"], flash["beta"]]

                    dew_isolated = self.phase_envelope_pt_mp(
                        z=z,
                        x_l0=x_l0,
                        w0=w0,
                        betas0=betas,
                        p0=flash["P"] - 5,
                        t0=flash["T"] + 5,
                        ns0=2 * len(z) + 3,
                        ds0=delta_bub_2ph,
                        max_points=max_points,
                    )

                    three_phase_envs.append((bubble_isolated, dew_isolated))

            stable_lines["2ph"] = {
                "dew": dew_line_stable,
                "bub": bub_line_stable,
                "liq": liq_line_stable,
            }

            return {
                "2ph": (dew_line, bub_line, liq_line),
                "DSP": dsps,
                "3ph": three_phase_envs,
                "2ph_stable": stable_lines["2ph"],
            }
        else:
            return {
                "2ph": {"dew": dew_line, "bub": bub_line},
            }

    # =========================================================================
    # Stability analysis
    # -------------------------------------------------------------------------
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
            - w: value of the test phase that minimizes the :math:`tm` function
            - tm: minimum value of the :math:`tm` function.
        dict
            All found minimum values of the :math:`tm` function and the
            corresponding test phase mole fractions.
            - w: all values of :math:`w` that minimize the :math:`tm` function
            - tm: all values found minima of the :math:`tm` function
        """
        (w_min, tm_min, all_mins) = yaeos_c.stability_zpt(
            id=self.id, z=z, p=pressure, t=temperature
        )

        all_mins_w = all_mins[:, : len(z)]
        all_mins = all_mins[:, -1]

        return {"w": w_min, "tm": tm_min}, {"tm": all_mins, "w": all_mins_w}

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

    # =========================================================================
    # Critical points and lines
    # -------------------------------------------------------------------------
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
        v0=0,
        t0=0,
        p0=0,
        stability_analysis=False,
        max_points=MAX_POINTS_ENVELOPES,
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
        v0: float, optional
            Initial guess for volume [L/mol]
        t0: float, optional
            Initial guess for temperature [K]
        p0: float, optional
            Initial guess for pressure [bar]
        max_points: int, optional
            Maximum number of points to calculate
        stop_pressure: float, optional
            Stop when reaching this pressure value
        """
        alphas, vs, ts, ps, *cep = yaeos_c.critical_line(
            self.id,
            ns=ns,
            ds0=ds0,
            a0=a0,
            v0=v0,
            t0=t0,
            p0=p0,
            s=s,
            stability_analysis=stability_analysis,
            z0=z0,
            zi=zi,
            max_points=max_points,
            stop_pressure=stop_pressure,
        )

        msk = ~np.isnan(ts)

        if stability_analysis:
            return {
                "a": alphas[msk],
                "T": ts[msk],
                "P": ps[msk],
                "V": vs[msk],
            }, {
                "x": cep[0],
                "y": cep[1],
                "P": cep[2],
                "Vx": cep[3],
                "Vy": cep[4],
                "T": cep[5],
            }
        else:
            return {
                "a": alphas[msk],
                "T": ts[msk],
                "P": ps[msk],
                "V": vs[msk],
            }

    def critical_line_liquid_liquid(
        self, z0=[0, 1], zi=[1, 0], pressure=2000, t0=500
    ):
        """Find the start of the Liquid-Liquid critical line of a binary.

        Parameters
        ----------
        z0: array_like
            Initial global mole fractions
        zi: array_like
            Final global mole fractions
        pressure: float
            Pressure [bar]
        t0: float
            Initial guess for temperature [K]
        """
        a, t, v = yaeos_c.find_llcl(
            self.id, z0=z0, zi=zi, p=pressure, tstart=t0
        )

        return a, t, v

    def __del__(self) -> None:
        """Delete the model from the available models list (Fortran side)."""
        yaeos_c.make_available_ar_models_list(self.id)
