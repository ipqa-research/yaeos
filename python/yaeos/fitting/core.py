"""Yaeos fitting core module."""

from typing import Callable

import numpy as np

import pandas as pd

from scipy.optimize import minimize

from yaeos.core import ArModel
from yaeos.fitting.solvers import solve_pt


KINDS = {
    "zT1": ["bubbleP", "bubble"],
    "zP1": ["bubbleT"],
    "zT0": ["dewP", "dew"],
    "zP0": ["dewT"],
    "PT": ["PT", "liquid-liquid"],
    "CP": ["critical", "CP"],
}


# =============================================================================
# Default error functions
# -----------------------------------------------------------------------------
def _pressure_error(Pexp, Pmodel):
    return (Pexp - Pmodel) ** 2 / Pexp


def _temperature_error(Texp, Tmodel):
    return (Texp - Tmodel) ** 2 / Texp


def _composition_error(zexp, zmodel):
    return np.abs(np.log(zmodel[0] / zexp[0])) + np.abs(
        np.log(zmodel[1] / zexp[1])
    )


class BinaryFitter:
    """BinaryFitter class.

    This class is used to fit binary interaction parameters to experimental
    data. The objective function is defined as the sum of the squared errors
    between the experimental data and the model predictions.

    Parameters
    ----------
    model_setter : Callable
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
        - x1: float, the mole fraction of component 1 (lightest component) in
        heavy phase.
        - y1: float, the mole fraction of component 1 (lightest component) in
        light phase.
        - T: float, the temperature in K
        - P: float, the pressure in bar
        - weight: float, optional, the weight of the data point (default is 1)
    verbose : bool, optional
        If True, print the objective function value and the optimization
    pressure_error : Callable, optional
        A function `f(Pexp, Pmodel)`that calculates the pressure error between 
        the experimental and model values. 
        The function should take the experimental and model
        values as arguments and return the error. If None, the default function
        is used.
    temperature_error : Callable, optional
        A function `f(Texp, Tmodel)`that calculates the temperature error
        between the experimental and model values. 
        The function should take the experimental and model
        values as arguments and return the error. If None, the default function
        is used.
    composition_error : Callable, optional
        A function `f(zexp, zmodel)`that calculates the composition error
        between the experimental and model values.  The function should take
        the experimental and model for mixture composition as arguments and
        return the error. If None, the default function is used.

    Attributes
    ----------
    get_model : Callable
        The function that returns the model object from the optimization
        parameters and the model_setter_args.
    get_model_args : tuple
        The arguments to pass to the model_setter function.
    data : pandas.DataFrame
        The experimental data.
    evaluations : dict
        A dictionary with the evaluations of the objective function. The keys
        are 'fobj', and 'x'. 'fobj' is the objective function value,
        'x' is the optimization parameters.
    """

    def __init__(
        self,
        model_setter: Callable,
        model_setter_args: tuple,
        data: pd.DataFrame,
        verbose: bool = False,
        pressure_error: Callable = None,
        temperature_error: Callable = None,
        composition_error: Callable = None,
    ) -> None:
        self.get_model = model_setter
        self.get_model_args = model_setter_args
        self.data = data
        self.verbose = verbose
        self.evaluations = {"fobj": [], "x": []}

        if pressure_error is not None:
            self.pressure_error = pressure_error
        else:
            self.pressure_error = _pressure_error

        if temperature_error is not None:
            self.temperature_error = temperature_error
        else:
            self.temperature_error = _temperature_error

        if composition_error is not None:
            self.composition_error = composition_error
        else:
            self.composition_error = _composition_error


    def objective_function(self, x_values) -> float:
        """
        Objective function to minimize when fitting interaction parameters.

        Parameters
        ----------
        x_values : array-like
            The interaction parameters to fit, 1D array-like.

        Returns
        -------
        float
            The value of the objective function, which is the sum of the
            squared relative errors between the experimental data and the model
            predictions.
        """
        model: ArModel = self.get_model(x_values, *self.get_model_args)
        data = self.data

        residuals = []
        objective_function_contributions = []

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
                if np.isnan(w) or w is None:
                    w = 1
            except KeyError:
                w = 1

            error_i = 0

            # =================================================================
            # Bubble point - Pressure
            # -----------------------------------------------------------------
            if row["kind"] in KINDS["zT1"]:
                sat = model.saturation_pressure(
                    x, kind="bubble", temperature=t, p0=p
                )
                error_i += self.pressure_error(p, sat["P"])
                residual = p - sat["P"]

            # =================================================================
            # Bubble point - Temperature
            # -----------------------------------------------------------------
            elif row["kind"] in KINDS["zP1"]:
                sat = model.saturation_temperature(
                    x, kind="bubble", pressure=p, t0=t
                )
                error_i += self.temperature_error(t, sat["T"])
                residual = t - sat["T"]

            # =================================================================
            # Dew point - Pressure
            # -----------------------------------------------------------------
            elif row["kind"] in KINDS["zT0"]:
                sat = model.saturation_pressure(
                    y, kind="dew", temperature=t, p0=p
                )
                error_i += self.pressure_error(p, sat["P"])
                residual = p - sat["P"]

            # =================================================================
            # Dew point - Temperature
            # -----------------------------------------------------------------
            elif row["kind"] in KINDS["zP0"]:
                sat = model.saturation_temperature(
                    y, kind="dew", pressure=p, t0=t
                )
                error_i += self.temperature_error(t, sat["T"])
                residual = t - sat["T"]

            # =================================================================
            # PT or Liquid-liquid equilibrium
            # -----------------------------------------------------------------
            elif row["kind"] in KINDS["PT"]:
                x1, y1 = solve_pt(model, row["P"], row["T"], row["kind"])

                residual = []
                if np.isnan(x[0]):
                    error_i += self.composition_error(y, [y1, 1 - y1])
                    residual.append(y[0] - y1)

                elif np.isnan(y[0]):
                    error_i += self.composition_error(x, [x1, 1 - x1])
                    residual.append(x[0] - x1)

                else:
                    error_i += self.composition_error(
                        x, [x1, 1 - x1]
                    ) + self.composition_error(y, [y1, 1 - y1])
                    residual.append(x[0] - x1)
                    residual.append(y[0] - y1)

            # =================================================================
            # Critical point error is calculated by finding the nearest
            # critical point in the critical line to the given critical
            # point in the data.
            # -----------------------------------------------------------------
            elif row["kind"] in KINDS["CP"]:
                cp = row
                distances = (cp["T"] - cl["T"]) ** 2 + (
                    (cp["P"] - cl["P"]) ** 2
                )

                if not np.isnan(cp["x1"]):
                    distances += (cp["x1"] - cl["a"]) ** 2
                nearest = np.argmin(distances)
                t_cl, p_cl, x1 = (
                    cl["T"][nearest],
                    cl["P"][nearest],
                    cl["a"][nearest],
                )
                error_i += self.temperature_error(cp["T"], t_cl)
                error_i += self.pressure_error(cp["P"], p_cl)
                error_i += self.composition_error(
                    [cp["x1"], 1 - cp["x1"]], [x1, 1 - x1]
                )
                residual = [
                    cp["T"] - t_cl,
                    cp["P"] - p_cl,
                    cp["x1"] - x1,
                ]

            else:
                raise ValueError(f"{row['kind']} is not a valid data kind.")

            if np.isnan(error_i) or np.isinf(error_i):
                # TODO make more robust PT solver to avoid infs
                error_i = row["P"]

            objective_function_contributions.append(error_i / len(data))
            residuals.append(residual)

            err += error_i * w

        # =====================================================================
        # Normalize the error and save the valuation
        # ---------------------------------------------------------------------
        err = err / len(data)

        self.evaluations["fobj"].append(err)
        self.evaluations["x"].append(x_values)
        self.evaluations["residuals"] = residuals
        self.evaluations["contributions"] = objective_function_contributions

        if self.verbose:
            print(len(self.evaluations["fobj"]), err, x_values)
        return err

    def fit(
        self, x0, bounds=None, method="Nelder-Mead", optimizer_options=None
        ):
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
            self.objective_function, x0=x0, bounds=bounds, method=method,
            options=optimizer_options
        )
        self._solution = sol

    @property
    def solution(self):
        """Return the optimization solution."""
        return self._solution

    @property
    def model(self):
        """Return the model with the fitted parameters."""
        if not hasattr(self, "_solution"):
            raise ValueError(
                "You must fit the model before getting the model."
            )
        return self.get_model(self._solution.x, *self.get_model_args)
