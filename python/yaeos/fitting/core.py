"""
"""

import numpy as np
from scipy.optimize import minimize

from yaeos.fitting.solvers import solve_PT


class BinaryFitter:

    def __init__(self, model_setter, model_setter_args, data, verbose=False):
        self._get_model = model_setter
        self._get_model_args = model_setter_args
        self.data = data
        self.verbose = verbose
        self.evaluations = {"fobj": [], "x": []}

    def objective_function(self, x_values):
        model = self._get_model(x_values, *self._get_model_args)
        data = self.data

        # ==============================================================
        # Calculate the critical line and remove the NaN values
        # --------------------------------------------------------------
        cl = model.critical_line(z0=[0, 1], zi=[1, 0], a0=1e-2, S=1e-2)
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
            w = row["weight"]
            error_i = 0

            # =================================================================
            # Bubble point
            # -----------------------------------------------------------------
            if row["kind"] == "bubble":
                sat = model.saturation_pressure(
                    x, kind="bubble", temperature=t, p0=p
                )
                error_i += (sat["P"] - p) ** 2 / p

                if not np.isnan(row["y1"]):
                    error_i += ((sat["y"][0] - y[0]) / y[0]) ** 2

            # =================================================================
            # Dew points
            # -----------------------------------------------------------------
            if row["kind"] == "dew":
                sat = model.saturation_pressure(
                    x, kind="dew", temperature=t, p0=p
                )
                error_i += (sat["P"] - p) ** 2 / p + (sat["y"][0] - y[0]) ** 2

            if row["kind"] == "PT" or row["kind"] == "liquid-liquid":
                x1, y1 = solve_PT(model, row["P"], row["T"])

                if np.isnan(x[0]):
                    error_i += (y1 - y[0]) ** 2

                elif np.isnan(y[0]):
                    error_i += (x1 - x[0]) ** 2

                else:
                    error_i += (x1 - row["x1"]) ** 2 + (y1 - row["y1"]) ** 2

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
                error_i += distances[nearest]

            err += error_i * w
        err = err / len(data)

        self.evaluations["fobj"].append(err)
        self.evaluations["x"].append(x_values)

        if self.verbose:
            print(err, x_values)
        return err

    def fit(self, x0, bounds, method="Nelder-Mead"):
        sol = minimize(
            self.objective_function, x0=x0, bounds=bounds, method=method
        )
        self._solution = sol
