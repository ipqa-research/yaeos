"""Global Phase Equilibria Calculations.

Module that implements the GPEC algorithm for calculation of GPEDs and its
derivatives to obtain isopleths, isotherms and isobars.
"""

import matplotlib.pyplot as plt

import numpy as np

from yaeos.core import ArModel

MAX_POINTS = 1000


class GPEC:
    """Global Phase Equilibria Calculation.

    This class implements the GPEC algorithm for calculation of global phase
    equilibria diagrams (GPEDs) and their usage to obtain isotherms, isobars
    and isopleths. It is designed to work with a given thermodynamic model
    that implements the necessary methods for phase equilibria calculations.

    Parameters
    ----------
    model : ArModel
        The thermodynamic model to be used for phase equilibria calculations.
    max_pressure : float, optional
        The maximum pressure for the calculations (default is 2500).
    max_points : int, optional
        The maximum number of points to be calculated in the critical line
        (default is 10000).
    stability_analysis : bool, optional
        Whether to perform stability analysis (default is True).
    step_21 : float, optional
        Step size for the critical line that starts from the almost pure second
        component (default is 1e-2).
    step_12 : float, optional
        Step size for the critical line that starts from the almost pure first
        component (default is 1e-5).
    x20 : float, optional
        Initial mole fraction of the second component for the critical line
        starting from the almost pure second component (default is 0.9999).
    x10 : float, optional
        Initial mole fraction of the first component for the critical line
        starting from the almost pure first component (default is 0.99999).

    Attributes
    ----------
    _z0 : np.ndarray
        Initial composition of the system, representing the almost pure second
        component.
    _zi : np.ndarray
        Initial composition of the system, representing the almost pure first
        component.
    _model : ArModel
        The thermodynamic model used for phase equilibria calculations.
    _pures : list
        List of pure saturation pressures for each component in the system.
    _cl21 : dict
        Critical line data starting from the almost pure second component.
    _cep21 : dict
        Critical endpoint data starting from the almost pure second component.
    _cl12 : dict or None
        Critical line data starting from the almost pure first component,
        or None if not applicable.
    _cep12 : dict or None
        Critical endpoint data starting from the almost pure first component,
        or None if not applicable.
    _cl_ll : dict
        Critical line data for liquid-liquid critical line.
    _cep_ll : dict
        Critical endpoint data for liquid-liquid critical line.

    Methods
    -------
    plot_gped()
        Plots the global phase equilibria diagram (GPED) based on the critical
        lines and pure saturation pressures.
    plot_pxy(temperature, a0=1e-5)
        Plots a Pxy phase diagram for the system at a given temperature.
    plot_txy(pressure, a0=1e-5)
        Plots a Txy phase diagram for the system at a given pressure.
    """

    def __init__(
        self,
        model: ArModel,
        max_pressure=2500,
        max_points=10000,
        stability_analysis=True,
        step_21=1e-2,
        step_12=1e-5,
        x20=0.9999,
        x10=0.99999,
    ):
        self._z0 = np.array([0, 1])
        self._zi = np.array([1, 0])
        self._model = model

        # Calculate the pure saturation pressures of each component and
        # save it as an internal variable
        psats = [model.pure_saturation_pressures(i) for i in [1, 2]]
        self._pures = psats

        # Calculate the critical line starting from the almost pure second
        # component.
        cl, cep = model.critical_line(
            z0=self._z0,
            zi=self._zi,
            ns=1,
            s=1 - x20,
            a0=1 - x20,
            ds0=step_21,
            stop_pressure=max_pressure,
            max_points=max_points,
            stability_analysis=stability_analysis,
        )

        self._cl21 = cl
        if np.isnan(cep["T"]):
            self._cep21 = None
        else:
            self._cep21 = cep

        # Check if the critical line did not reach to the pure first component.
        # if not, calculate the critical line starting from the almost pure
        # first component. It is important to make small steps because this
        # kind of line can be pretty short.
        if not np.isnan(cep["T"]) or (
            abs(cl["T"][-1] - psats[0]["T"][-1]) > 10
            and abs(cl["P"][-1] - psats[0]["P"][-1] > 10)
        ):
            cl, cep = model.critical_line(
                z0=self._z0,
                zi=self._zi,
                ns=1,
                s=x10,
                a0=x10,
                ds0=-step_12,
                stop_pressure=max_pressure,
                max_points=max_points,
                stability_analysis=stability_analysis,
            )

            self._cl12 = cl

            if np.isnan(cep["T"]):
                self._cep12 = None
            else:
                self._cep12 = cep
        else:
            self._cl12 = None
            self._cep12 = None

        a, temp, vol = model.critical_line_liquid_liquid(
            z0=self._z0, zi=self._zi, pressure=max_pressure, t0=500
        )

        cl, cep = model.critical_line(
            z0=self._z0,
            zi=self._zi,
            ns=4,
            s=np.log(max_pressure),
            a0=a,
            v0=vol,
            t0=temp,
            p0=max_pressure,
            ds0=-1e-1,
            stop_pressure=max_pressure * 1.1,
            max_points=max_points,
            stability_analysis=stability_analysis,
        )

        if len(cl["a"]) > 0:
            self._cl_ll = cl
        else:
            self._cl_ll = None

        if np.isnan(cep["T"]):
            self._cep_ll = None
        else:
            self._cep_ll = cep

        # Calculate the three-phase lines from the critical endpoints
        if self._cep12:
            self._llv_12 = model.binary_llv_from_cep(self._cep12)
        if self._cep21:
            self._llv_21 = model.binary_llv_from_cep(self._cep21)
        if self._cep_ll:
            self._llv_ll = model.binary_llv_from_cep(self._cep_ll)

    def plot_gped(self):
        """Plot the global phase equilibria diagram (GPED).

        This method plots the critical lines and pure saturation pressures
        for the components in the system. It visualizes the phase behavior of
        the system across different temperatures and pressures.
        """
        for pure in self._pures:
            plt.plot(pure["T"], pure["P"], color="green")

        plt.plot(self._cl21["T"], self._cl21["P"], color="black")
        if self._cl12:
            plt.plot(self._cl12["T"], self._cl12["P"], color="black")
        if self._cl_ll:
            plt.plot(self._cl_ll["T"], self._cl_ll["P"], color="black")
        if self._cep12:
            plt.plot(self._llv_12["T"], self._llv_12["P"], color="purple", linestyle="--")
        if self._cep21:
            plt.plot(self._llv_21["T"], self._llv_21["P"], color="purple", linestyle="--")
        if self._cep_ll:
            plt.plot(self._llv_ll["T"], self._llv_ll["P"], color="purple", linestyle="--")

        plt.plot([], [], color="green", label="Pure saturation pressure")
        plt.plot([], [], color="black", label="Critical line")
        plt.plot([], [], color="purple", linestyle="--", label="Three-phase line")
        plt.legend(frameon=False)
        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (bar)")

    def calc_pxy(
        self,
        temperature,
        x10=0.9999,
        x20=0.9999,
        dx10=1e-3,
        dx20=-1e-3,
        dxll0=1e-4,
    ):
        """Calculate a Pxy phase diagram.

        This method calculates the Pxy phase diagram for the system at a given
        temperature. It uses the critical lines and pure saturation pressures
        to determine the phase behavior of the system.

        Parameters
        ----------
        temperature : float
            The temperature at which to calculate the Pxy phase diagram.
        x10 : float, optional
            Initial molar fraction of the first component for the line that
            starts from the almost pure first component (default is 1e-5).
        x20 : float, optional
            Initial molar fraction of the second component for the line that
            starts from the almost pure second component (default is 1e-5).
        dx10 : float, optional
            Step size for the line that starts from the almost pure first
            component (default is 1e-3).
        dx20 : float, optional
            Step size for the line that starts from the almost pure second
            component (default is 1e-3).
        dxll0 : float, optional
            Step size for the line that starts from the liquid-liquid critical
            point (default is 1e-4).

        Returns
        -------
        list
            A list containing the calculated Pxy phase diagrams.
            - Starting from the almost pure second component.
            - Starting from the almost pure first component. If below its Tc.
            - Starting from liquid-liquid critical line, if it exists.
        """
        psat_1, psat_2 = self._pures
        critical_temperatures = max(psat_1["T"]), max(psat_2["T"])

        px_12 = px_21 = px_ll = None

        # Below saturation temperature of light component
        loc = np.argmin(abs(psat_2["T"] - temperature))
        p0 = psat_2["P"][loc]
        px_21 = self._model.phase_envelope_px(
            self._z0,
            self._zi,
            temperature=temperature,
            kind="bubble",
            p0=p0,
            a0=x10,
            ds0=dx10,
            max_points=MAX_POINTS,
        )

        # Below saturation temperature of heavy component
        if temperature < critical_temperatures[1]:
            loc = np.argmin(abs(psat_1["T"] - temperature))
            p0 = psat_1["P"][loc]
            px_12 = self._model.phase_envelope_px(
                self._z0,
                self._zi,
                temperature=temperature,
                kind="bubble",
                p0=p0,
                a0=1 - x20,
                ds0=dx20,
                max_points=MAX_POINTS,
            )

        if self._cl_ll:
            loc = np.argmin(abs(self._cl_ll["T"] - temperature))
            p0, t = self._cl_ll["P"][loc], self._cl_ll["T"][loc]
            if abs(t - temperature) < 1 or temperature > self._cep_ll["T"]:

                a = self._cl_ll["a"][loc]
                z = a * self._zi + (1 - a) * self._z0
                x_l0 = [z.copy()]

                x_l0[0][0] += 1e-5
                x_l0[0][1] -= 1e-5
                w0 = z.copy()
                w0[0] -= 1e-5
                w0[1] += 1e-5

                px_ll = self._model.phase_envelope_px_mp(
                    z0=self._z0,
                    zi=self._zi,
                    t=temperature,
                    kinds_x=["liquid"],
                    kind_w="liquid",
                    x_l0=x_l0,
                    w0=w0,
                    betas0=[1],
                    p0=p0,
                    alpha0=a + 1e-5,
                    ns0=len(z) + 3,
                    ds0=dxll0,
                    max_points=MAX_POINTS,
                )

        pxs = [px_12, px_21, px_ll]

        return pxs

    def plot_pxy(
        self,
        temperature,
        x10=1e-5,
        x20=1e-5,
        dx10=1e-3,
        dx20=1e-3,
        dxll0=1e-4,
        **plot_kwargs,
    ):
        """Plot a Pxy phase diagram.

        This method plots the Pxy phase diagram for the system at a given
        temperature. It visualizes the phase behavior of the system across
        different compositions of the components.

        Parameters
        ----------
        temperature : float
            The temperature at which to plot the Pxy phase diagram.
        x10 : float, optional
            Initial molar fraction of the first component for the line that
            starts from the almost pure first component (default is 1e-5).
        x20 : float, optional
            Initial molar fraction of the second component for the line that
            starts from the almost pure second component (default is 1e-5).
        dx10 : float, optional
            Step size for the line that starts from the almost pure first
            component (default is 1e-3).
        dx20 : float, optional
            Step size for the line that starts from the almost pure second
            component (default is 1e-3).
        dxll0 : float, optional
            Step size for the line that starts from the liquid-liquid critical
            point (default is 1e-4).
        plot_kwargs : dict, optional
            Additional keyword arguments to be passed to the plot function.
        """
        pxs = self.calc_pxy(
            temperature, x10=x10, x20=x20, dx10=dx10, dx20=dx20, dxll0=dxll0
        )

        for px in pxs:
            if px:
                plt.plot(
                    px.main_phases_compositions[:, 0, 0],
                    px["P"],
                    **plot_kwargs,
                )
                plt.plot(
                    px.reference_phase_compositions[:, 0],
                    px["P"],
                    **plot_kwargs,
                )

        plt.xlabel(r"$x_1$, $y_1$")
        plt.ylabel("Pressure (bar)")

    def calc_txy(
        self,
        pressure,
        y10=0.9999,
        y20=0.9999,
        dy10=-1e-3,
        dy20=1e-3,
        dyll0=1e-4,
    ):
        """Calculate a Txy phase diagram.

        This method calculates the Txy phase diagram for the system at a given
        pressure. It uses the critical lines and pure saturation pressures to
        determine the phase behavior of the system.

        Parameters
        ----------
        pressure : float
            The pressure at which to calculate the Txy phase diagram.
        y10 : float, optional
            Initial molar fraction of the first component for the line that
            starts from the almost pure first component (default is 1e-5).
        y20 : float, optional
            Initial molar fraction of the second component for the line that
            starts from the almost pure second component (default is 1e-5).
        dy10 : float, optional
            Step size for the line that starts from the almost pure first
            component (default is -1e-3).
        dy20 : float, optional
            Step size for the line that starts from the almost pure second
            component (default is 1e-3).
        dyll0 : float, optional
            Step size for the line that starts from the liquid-liquid critical
            point (default is 1e-4).

        Returns
        -------
        list
            A list containing the calculated Txy envelopes. They are sorted
            in the order of:
            - Starting from the almost pure first component.
            - Starting from the almost pure second component.
            - Starting from liquid-liquid critical line if it exists.
        """
        psat_1, psat_2 = self._pures

        tx_12 = tx_21 = tx_ll = None

        # Below critical pressure of heavy component
        if pressure < psat_2["P"][-1]:
            loc = np.argmin(abs(psat_2["P"] - pressure))
            t0 = psat_2["T"][loc]

            tx_21 = self._model.phase_envelope_tx(
                self._z0,
                self._zi,
                pressure=pressure,
                kind="dew",
                t0=t0,
                a0=1 - y20,
                ds0=1e-5,
                max_points=MAX_POINTS,
            )
        if pressure < psat_1["P"][-1]:
            # Below critical pressure of the light component
            loc = np.argmin(abs(psat_1["P"] - pressure))
            t0 = psat_1["T"][loc]

            tx_12 = self._model.phase_envelope_tx(
                self._z0,
                self._zi,
                pressure=pressure,
                kind="bubble",
                t0=t0,
                a0=y10,
                ds0=dy10,
                max_points=MAX_POINTS,
            )
        if self._cl_ll:
            loc = np.argmin(abs(self._cl_ll["P"] - pressure))
            t0, p = self._cl_ll["T"][loc], self._cl_ll["P"][loc]
            if abs(p - pressure) < 1 or pressure > psat_1["P"][-1]:

                a = self._cl_ll["a"][loc]
                z = a * self._zi + (1 - a) * self._z0
                x_l0 = [z.copy()]

                x_l0[0][0] += 1e-5
                x_l0[0][1] -= 1e-5
                w0 = z.copy()
                w0[0] -= 1e-5
                w0[1] += 1e-5

                tx_ll = self._model.phase_envelope_tx_mp(
                    z0=self._z0,
                    zi=self._zi,
                    p=pressure,
                    kinds_x=["liquid"],
                    kind_w="liquid",
                    x_l0=x_l0,
                    w0=w0,
                    betas0=[1],
                    t0=t0,
                    alpha0=a + 1e-5,
                    ns0=len(z) + 3,
                    ds0=dyll0,
                    max_points=MAX_POINTS,
                )

        txs = [tx_12, tx_21, tx_ll]

        return txs

    def plot_txy(
        self,
        pressure,
        y10=0.9999,
        y20=0.9999,
        dy10=-1e-3,
        dy20=1e-3,
        dyll0=1e-4,
        **plot_kwargs,
    ):
        """Plot a Txy phase diagram.

        This method plots the Txy phase diagram for the system at a given
        pressure.

        Parameters
        ----------
        pressure : float
            The pressure at which to plot the Txy phase diagram.
        y10 : float, optional
            Initial molar fraction of the first component for the line that
            starts from the almost pure first component (default is 0.9999).
        y20 : float, optional
            Initial molar fraction of the second component for the line that
            starts from the almost pure second component (default is 0.9999).
        dy10 : float, optional
            Step size for the line that starts from the almost pure first
            component (default is -1e-3).
        dy20 : float, optional
            Step size for the line that starts from the almost pure second
            component (default is 1e-3).
        dyll0 : float, optional
            Step size for the line that starts from the liquid-liquid critical
            point (default is 1e-4).
        plot_kwargs : dict, optional
            Additional keyword arguments to be passed to the plot function.
        """
        txs = self.calc_txy(
            pressure, y10=y10, y20=y20, dy10=dy10, dy20=dy20, dyll0=dyll0
        )

        for tx in txs:
            if tx:
                plt.plot(
                    tx.main_phases_compositions[:, 0, 0],
                    tx["T"],
                    **plot_kwargs,
                )
                plt.plot(
                    tx.reference_phase_compositions[:, 0],
                    tx["T"],
                    **plot_kwargs,
                )

        plt.xlabel(r"$x_1$, $y_1$")
        plt.ylabel("Temperature (K)")
