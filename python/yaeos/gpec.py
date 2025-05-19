"""Global Phase Equilibria Calculations.

Module that implements the GPEC algorithm for calculation of GPEDs and its
derivatives to obtain isopleths, isotherms and isobars.
"""

import matplotlib.pyplot as plt

import numpy as np

from yaeos.core import ArModel


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
        self._z0 = [0, 1]
        self._zi = [1, 0]
        self._model = model

        psats = [model.pure_saturation_pressures(i) for i in [1, 2]]

        self._pures = psats

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

        if (
            not np.isnan(cep["T"])
            or (
                abs(cl["T"][-1] - psats[0]["T"][-1]) > 10
                and abs(cl["P"][-1] - psats[0]["P"][-1] > 10)
            )
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
            plt.plot(self._cl12["T"], self.cl12["P"], color="black")
        if self._cl_ll:
            plt.plot(self._cl_ll["T"], self._cl_ll["P"], color="black")