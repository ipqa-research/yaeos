"""Global Phase Equilibria Calculations.

Module that implements the GPEC algorithm for calculation of GPEDs and its
derivatives to obtain isopleths, isotherms and isobars.
"""

import matplotlib.pyplot as plt

import numpy as np

from yaeos.core import ArModel

MAX_POINTS = 10000


class GPEC:

    def __init__(
        self,
        model: ArModel,
        max_pressure=2500,
        max_points=10000,
        stability_analysis=True,
        step_21=1e-2,
        step_12=1e-5,
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
        diff = 1e-3
        cl, cep = model.critical_line(
            z0=self._z0,
            zi=self._zi,
            ns=1,
            s=diff,
            a0=diff,
            ds0=step_21,
            stop_pressure=max_pressure,
            max_points=max_points,
            stability_analysis=stability_analysis,
        )

        self._cl21 = cl
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
                s=1 - diff / 10,
                a0=1 - diff / 10,
                ds0=-step_12,
                stop_pressure=max_pressure,
                max_points=max_points,
                stability_analysis=stability_analysis,
            )

            self._cl12 = cl
            self._cep12 = cep
        else:
            self._cl12 = None
            self._cep12 = None

        a, T, V = model.critical_line_liquid_liquid(
            z0=self._z0, zi=self._zi, pressure=max_pressure, t0=500
        )

        cl, cep = model.critical_line(
            z0=self._z0,
            zi=self._zi,
            ns=4,
            s=np.log(max_pressure),
            a0=a,
            v0=V,
            t0=T,
            p0=max_pressure,
            ds0=-1e-1,
            stop_pressure=max_pressure * 1.1,
            max_points=max_points,
            stability_analysis=stability_analysis,
        )

        self._cl_ll = cl
        self._cep_ll = cep

    def plot_gped(self):
        for pure in self._pures:
            plt.plot(pure["T"], pure["P"], color="green")

        plt.plot(self._cl21["T"], self._cl21["P"], color="black")
        if self._cl12:
            plt.plot(self._cl12["T"], self._cl12["P"], color="black")
        if self._cl_ll:
            plt.plot(self._cl_ll["T"], self._cl_ll["P"], color="black")

        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (bar)")

    def plot_pxy(self, temperature, a0=1e-5):
        """Plot a Pxy phase diagram"""
        psat_1, psat_2 = self._pures

        px_12 = px_21 = px_iso = None

        # Below saturation temperature of light component
        loc = np.argmin(abs(psat_2["T"] - temperature))
        p0 = psat_2["P"][loc]
        px_21 = self._model.phase_envelope_px(
            self._z0,
            self._zi,
            temperature=temperature,
            kind="bubble",
            p0=p0,
            a0=a0,
            ds0=1e-3,
            max_points=MAX_POINTS,
        )

        pxs = [px_12, px_21, px_iso]

        for px in pxs:
            if px:
                plt.plot(px.main_phases_compositions[:, 0, 0], px["P"])
                plt.plot(px.reference_phase_compositions[:, 0], px["P"])

        plt.xlabel(r"$x_1$, $y_1$")
        plt.ylabel("Pressure (bar)")

        return pxs

    def plot_txy(self, pressure, a0=1e-5):
        """Plot a Txy phase diagram"""
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
                a0=a0,
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
                a0=1 - a0,
                ds0=-1e-5,
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
                    ds0=1e-4,
                    max_points=MAX_POINTS,
                )

        txs = [tx_12, tx_21, tx_ll]

        for tx in txs:
            if tx:
                plt.plot(tx.main_phases_compositions[:, 0, 0], tx["T"])
                plt.plot(tx.reference_phase_compositions[:, 0], tx["T"])

        plt.xlabel(r"$x_1$, $y_1$")
        plt.ylabel("Temperature (K)")

        return txs
