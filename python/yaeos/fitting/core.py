"""
"""

import numpy as np
from scipy.optimize import minimize

from .solvers import solve_PT


class BinaryFitter:

    def __init__(
            self, model_setter,
            model_setter_args,
            data, verbose=False):
        self._get_model = model_setter
        self._get_model_args = model_setter_args
        self.data = data
        self.verbose = verbose

        self.critical_points = len(data[data["kind"] == "critical"])
        self.bubble_points = len(data[data["kind"] == "bubble"])

        self._critical_weight = (len(data) - self.critical_points) / len(data)
        self._bubble_weigth = (len(data) - self.bubble_points) / len(data)

    def objective_function(self, x_values):
        model = self._get_model(x_values, *self._get_model_args)
        data = self.data

        w_crit = self._critical_weight
        w_bub = self._bubble_weigth

        # ==============================================================
        # Calculate the critical line and remove the NaN values
        # --------------------------------------------------------------
        cl = model.critical_line(z0=[0, 1], zi=[1, 0])
        msk = ~np.isnan(cl["T"])
        cl["T"] = cl["T"][msk]
        cl["P"] = cl["P"][msk]
        cl["a"] = cl["a"][msk]

        err = 0

        for _, row in data.iterrows():
            x = [row["x1"], 1 - row["x1"]]
            y = [row["y1"], 1 - row["y1"]]
            t = row["T"]
            p = row["P"]

            # =================================================================
            # Bubble point
            # -----------------------------------------------------------------
            if row["kind"] == "bubble":
                sat = model.saturation_pressure(
                    x, kind="bubble", temperature=t, p0=p
                )
                err += (sat["P"] - p) ** 2 / p + (sat["y"][0] - y[0]) ** 2
                err *= w_bub

            if row["kind"] == "PT":
                x1, y1 = solve_PT(model, row["P"], row["T"])
                err += (x1 - row["x1"]) ** 2 + (y1 - row["y1"]) ** 2

            # =================================================================
            # Critical point error is calculated by finding the nearest
            # critical point in the critical line to the given critical
            # point in the data.
            # -----------------------------------------------------------------
            if row["kind"] == "critical":
                cp = row
                distances = (
                    (cp["T"] - cl["T"]) ** 2
                    + ((cp["P"] - cl["P"]) ** 2)
                    + ((cp["x1"] - cl["a"]) ** 2)
                )
                nearest = np.argmin(distances)
                err += distances[nearest] * w_crit

        err = err/len(data)
        if self.verbose:
            print(err)
        return err

    def fit(self, x0, bounds):
        sol = minimize(
            self.objective_function,
            x0=x0,
            bounds=bounds,
            method="Nelder-Mead"
        )
        self._solution = sol
