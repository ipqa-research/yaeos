"""Yaeos fitting core module."""

import numpy as np

from scipy.optimize import minimize

from yaeos.fitting.solvers import solve_pt


class BinaryFitter:
    """BinaryFitter class.

    This class is used to fit binary interaction parameters to experimental
    data. The objective function is defined as the sum of the squared errors
    between the experimental data and the model predictions.

    Parameters
    ----------
    model_setter : callable
        A function that returns a model object. The function should take the
        optimization parameters as the first argument and any other arguments
        as the following arguments.
    model_setter_args : tuple
        A tuple with the arguments to pass to the model_setter function.
    data : pandas.DataFrame
        A DataFrame with the experimental data.
        The DataFrame should have the following columns:
        - kind: str, the kind of data point (bubble, dew, liquid-liquid, PT,
        critical)
        - x1: float, the mole fraction of component 1
        - y1: float, the mole fraction of component 1
        - T: float, the temperature in K
        - P: float, the pressure in bar
    verbose : bool, optional
        If True, print the objective function value and the optimization
    """

    def __init__(self, model_setter, model_setter_args, data, verbose=False):
        self._get_model = model_setter
        self._get_model_args = model_setter_args
        self.data = data
        self.verbose = verbose
        self.evaluations = {"fobj": [], "x": [], "cl": []}

    def objective_function(self, x_values):
        """
        Objective function to minimize when fitting interaction parameters.

        Parameters
        ----------
        x_values : array-like
            The interaction parameters to fit.
        """

        def pressure_error(Pexp, Pmodel):
            return (Pexp - Pmodel) ** 2 / Pexp

        def composition_error(zexp, zmodel):
            return np.abs(np.log(zmodel[0] / zexp[0])) + np.abs(
                np.log(zmodel[1] / zexp[1])
            )

        def temperature_error(Texp, Tmodel):
            return (Texp - Tmodel) ** 2 / Texp

        model = self._get_model(x_values, *self._get_model_args)
        data = self.data

        # =====================================================================
        # Calculate the critical line starting from the heavy component
        # ---------------------------------------------------------------------
        cp_msk = data["kind"] == "critical"
        if len(data[cp_msk]) > 0:
            cl = model.critical_line(
                z0=[0, 1],
                zi=[1, 0],
                a0=1e-2,
                s=1e-2,
                ds0=1e-3,
                max_points=15000,
            )

        err = 0

        for _, row in data.iterrows():
            x = [row["x1"], 1 - row["x1"]]
            y = [row["y1"], 1 - row["y1"]]
            t = row["T"]
            p = row["P"]

            try:
                w = row["weight"]
            except KeyError:
                w = 1

            error_i = 0

            # =================================================================
            # Bubble point
            # -----------------------------------------------------------------
            if row["kind"] == "bubble":
                sat = model.saturation_pressure(
                    x, kind="bubble", temperature=t, p0=p
                )
                error_i += pressure_error(p, sat["P"])

            # =================================================================
            # Dew points
            # -----------------------------------------------------------------
            if row["kind"] == "dew":
                sat = model.saturation_pressure(
                    x, kind="dew", temperature=t, p0=p
                )

                error_i += pressure_error(p, sat["P"])

            if row["kind"] == "PT" or row["kind"] == "liquid-liquid":
                x1, y1 = solve_pt(model, row["P"], row["T"], row["kind"])

                if np.isnan(x[0]):
                    error_i += composition_error(y, [y1, 1 - y1])

                elif np.isnan(y[0]):
                    error_i += composition_error(x, [x1, 1 - x1])

                else:
                    error_i += composition_error(
                        x, [x1, 1 - x1]
                    ) + composition_error(y, [y1, 1 - y1])

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
                t_cl, p_cl, x1 = (
                    cl["T"][nearest],
                    cl["P"][nearest],
                    cl["a"][nearest],
                )
                error_i += temperature_error(cp["T"], t_cl)
                error_i += pressure_error(cp["P"], p_cl)
                error_i += composition_error(
                    [cp["x1"], 1 - cp["x1"]], [x1, 1 - x1]
                )

            if np.isnan(error_i):
                error_i = row["P"]
            err += error_i * w

        # =====================================================================
        # Normalize the error and save the valuation
        # ---------------------------------------------------------------------
        err = err / len(data)

        self.evaluations["fobj"].append(err)
        self.evaluations["x"].append(x_values)

        if self.verbose:
            print(err, x_values)
        return err

    def fit(self, x0, bounds, method="Nelder-Mead"):
        """Fit the model to the data.

        Fit the model to the data using the objective function defined in
        the objective_function method. The optimization is performed using
        the scipy.optimize.minimize function.
        The optimization result is stored in the `.solution` property. Which

        Parameters
        ----------
        x0 : array-like
            Initial guess for the fitting parameters.
        bounds : array-like
            Bounds for the fitting parameters.
        method : str, optional
            The optimization method to use. Default is 'Nelder-Mead'.

        Returns
        -------
        None
        """
        sol = minimize(
            self.objective_function, x0=x0, bounds=bounds, method=method
        )
        self._solution = sol

    @property
    def solution(self):
        return self._solution
