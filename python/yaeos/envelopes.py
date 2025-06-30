"""Envelopes.

This module contains the classes that wrap the data structures used to
represent different kinds of phase envelopes.
"""

from IPython.display import display

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd


class PTEnvelope:
    """PTEnvelope.
    This class represents a pressure-temperature envelope.
    Attributes
    ----------
    global_composition : np.ndarray
        The global composition of the system. Shape is (n_components,).
    main_phases_compositions : np.ndarray
        The compositions of the main phases along the envelope.
        Shape is (n_points, n_phases, n_components).
    reference_phase_compositions : np.ndarray
        The compositions of the reference phase along the envelope.
        Shape is (n_points, n_components).
    main_phases_molar_fractions : np.ndarray
        The molar fractions of the main phases along the envelope.
        Shape is (n_points, n_phases).
    pressures : np.ndarray
        The pressures along the envelope. [bar]
    temperatures : np.ndarray
        The temperatures along the envelope. [K]
    iterations : np.ndarray
        The number of iterations taken to compute the envelope at each point.
        Shape is (n_points,).
    specified_variable : np.ndarray
        The specified variable used to compute the envelope at each point.
        Shape is (n_points,).
    reference_phase_kinds : np.ndarray
        The kinds of the reference phase at each point.
        Shape is (n_points,).
    main_phases_kinds : np.ndarray
        The kinds of the main phases at each point.
        Shape is (n_points, n_phases).
    cp : list
        A list of lists containing the indices of the critical points for each
        phase. Each sublist corresponds to a phase and contains the indices of
        the critical points in the `temperatures` and `pressures` arrays.
    df : pd.DataFrame
        A DataFrame containing the data of the envelope. The columns are:
        - 'T': Temperatures along the envelope.
        - 'P': Pressures along the envelope.
        - 'x_i^j': Compositions of the main phases, where `i` is the component index and `j` is the phase index.
        - 'w_i': Compositions of the reference phase, where `i` is the component index.
        - 'beta^j': Molar fractions of the main phases, where `j` is the phase index.
    """

    def __init__(
        self,
        global_composition,
        main_phases_compositions,
        reference_phase_compositions,
        reference_phase_kinds,
        main_phases_kinds,
        main_phases_molar_fractions,
        pressures,
        temperatures,
        iterations,
        specified_variable,
        critical_pressures,
        critical_temperatures,
    ):

        msk = ~np.isnan(pressures)
        msk_cp = ~np.isnan(critical_pressures)
        self.number_of_components = len(global_composition)
        self.number_of_phases = main_phases_compositions.shape[1]
        self.global_composition = global_composition
        self.main_phases_compositions = main_phases_compositions[msk, :, :]
        self.reference_phase_compositions = reference_phase_compositions[
            msk, :
        ]
        self.main_phases_molar_fractions = main_phases_molar_fractions[msk]
        self.pressures = pressures[msk]
        self.temperatures = temperatures[msk]
        self.iterations = iterations[msk]
        self.specified_variable = specified_variable[msk]
        self.reference_phase_kinds = reference_phase_kinds[msk]
        self.main_phases_kinds = main_phases_kinds[msk]
        self.critical_pressures = critical_pressures[msk_cp]
        self.critical_temperatures = critical_temperatures[msk_cp]

        df = pd.DataFrame()

        df["T"] = self.temperatures
        df["P"] = self.pressures

        for i in range(self.number_of_components):
            for j in range(self.number_of_phases):
                df[f"x_{i+1}^{j+1}"] = self.main_phases_compositions[:, j, i]
            df[f"w_{i+1}"] = self.reference_phase_compositions[:, i]

        for i in range(self.number_of_phases):
            df[f"beta^{i+1}"] = self.main_phases_molar_fractions[:, i]

        self.df = df

    def plot(self, **plot_kwargs):
        if "ax" in plot_kwargs:
            ax = plot_kwargs["ax"]
            del plot_kwargs["ax"]
        else:
            ax = plt.gca()
        ax.plot(self.temperatures, self.pressures, **plot_kwargs)
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("Pressure [bar]")
        ax.scatter(
            self.critical_temperatures, self.critical_pressures, color="black"
        )

    def __getitem__(self, key):
        if "key" in self.__dict__:
            return self.__dict__["key"]
        elif isinstance(key, np.ndarray) or isinstance(key, list):
            return PTEnvelope(
                global_composition=self.global_composition,
                main_phases_compositions=self.main_phases_compositions[key],
                reference_phase_compositions=self.reference_phase_compositions[
                    key
                ],
                main_phases_molar_fractions=self.main_phases_molar_fractions[
                    key
                ],
                main_phases_kinds=self.main_phases_kinds[key],
                reference_phase_kinds=self.reference_phase_kinds[key],
                pressures=self.pressures[key],
                temperatures=self.temperatures[key],
                iterations=self.iterations[key],
                specified_variable=self.specified_variable[key],
                critical_pressures=self.critical_pressures,
                critical_temperatures=self.critical_temperatures,
            )
        elif key == "T":
            return self.temperatures
        elif key == "Tc":
            return self.critical_temperatures
        elif key == "Pc":
            return self.critical_pressures
        elif key == "P":
            return self.pressures
        elif key == "z":
            return self.global_composition
        elif key == "x":
            return self.main_phases_compositions
        elif key == "w":
            return self.reference_phase_compositions

    def __repr__(self):
        display(self.df)
        return ""

    def __len__(self):
        return len(self.temperatures)


class PXEnvelope:
    """PXEnvelope.

    This class represents a pressure-composition envelope.

    Attributes
    ----------
    temperature : float
        The temperature at which the envelope is defined. [K]
    global_composition_0 : np.ndarray
        The global composition at the point where :math:`\alpha = 0`.
    global_composition_i : np.ndarray
        The global composition at the point where :math:`\alpha = 1`.
    main_phases_compositions : np.ndarray
        The compositions of the main phases along the envelope.
        Shape is (n_points, n_phases, n_components).
    reference_phase_compositions : np.ndarray
        The compositions of the reference phase along the envelope.
        Shape is (n_points, n_components).
    main_phases_molar_fractions : np.ndarray
        The molar fractions of the main phases along the envelope.
        Shape is (n_points, n_phases).
    pressures : np.ndarray
        The pressures along the envelope. [bar]
    alphas : np.ndarray
        The molar fraction of the `global_composition_i`, :math:`\alpha`.
        Shape is (n_points,).
    critical_pressures : np.ndarray
        The critical pressures of the envelope. Shape is (n_critical_points,).
    critical_alphas : np.ndarray
        The molar fractions of the `global_composition_i` at the critical
        points, :math:`\alpha`.  Shape is (n_critical_points,).
    iterations : np.ndarray
        The number of iterations taken to compute the envelope at each point.
        Shape is (n_points,).
    specified_variable : np.ndarray
        The specified variable used to compute the envelope at each point.
        Shape is (n_points,).
    df : pd.DataFrame
        A DataFrame containing the data of the envelope. The columns are:
        - 'alpha': Molar fraction of the `global_composition_i`.
        - 'P': Pressures along the envelope.
        - 'x_i^j': Compositions of the main phases, where `i`
        is the component index and `j` is the phase index.
        - 'w_i': Compositions of the reference phase, where `i` is the
        component index.
        - 'beta^j': Molar fractions of the main phases, where `j`
        is the phase index.
    """

    def __init__(
        self,
        temperature,
        global_composition_0,
        global_composition_i,
        main_phases_compositions,
        reference_phase_compositions,
        main_phases_molar_fractions,
        pressures,
        alphas,
        iterations,
        specified_variable,
        critical_pressures,
        critical_alphas,
        main_phases_kinds,
        reference_phase_kinds,
    ):

        msk = ~np.isnan(pressures)
        msk_cp = ~np.isnan(critical_pressures)

        self.temperature = temperature
        self.number_of_components = len(global_composition_0)
        self.number_of_phases = main_phases_compositions.shape[1]
        self.global_composition_0 = global_composition_0
        self.global_composition_i = global_composition_i
        self.main_phases_compositions = main_phases_compositions[msk, :, :]
        self.reference_phase_compositions = reference_phase_compositions[
            msk, :
        ]
        self.main_phases_molar_fractions = main_phases_molar_fractions[msk]
        self.pressures = pressures[msk]
        self.alphas = alphas[msk]
        self.iterations = iterations[msk]
        self.specified_variable = specified_variable[msk]
        self.main_phases_kinds = main_phases_kinds[msk]
        self.reference_phase_kinds = reference_phase_kinds[msk]
        self.critical_pressures = critical_pressures[msk_cp]
        self.critical_alphas = critical_alphas[msk_cp]

        df = pd.DataFrame()

        df["alpha"] = self.alphas
        df["P"] = self.pressures

        for i in range(self.number_of_components):
            for j in range(self.number_of_phases):
                df[f"x_{i+1}^{j+1}"] = self.main_phases_compositions[:, j, i]
            df[f"w_{i+1}"] = self.reference_phase_compositions[:, i]

        for i in range(self.number_of_phases):
            df[f"beta^{i+1}"] = self.main_phases_molar_fractions[:, i]

        self.df = df

    def plot(self, **plot_kwargs):
        """Plot the envelope.

        Plot the envelope in a matplotlib axis. Using the natural variables of
        the envelope (:math:`alpha` and pressure)  If the `ax` keyword argument
        is provided, it will be used as the axis to plot on. Otherwise, the
        current axis will be used.

        Parameters
        ----------
        plot_kwargs : dict
            Keyword arguments to pass to the `plot` method of the axis.
            This can include line style, color, etc.
        """
        if "ax" in plot_kwargs:
            ax = plot_kwargs["ax"]
            del plot_kwargs["ax"]
        else:
            ax = plt.gca()
        ax.plot(self.alphas, self.alphas, **plot_kwargs)
        ax.set_xlabel(r"$\alpha$")
        ax.set_ylabel("Pressure [bar]")
        for cp in self.cp:
            ax.scatter(self.pressures[cp], self.alphas[cp], color="black")

    def __getitem__(self, key):
        if "key" in self.__dict__:
            return self.__dict__["key"]
        elif isinstance(key, np.ndarray) or isinstance(key, list):
            return PXEnvelope(
                global_composition_0=self.global_composition_0,
                global_composition_i=self.global_composition_i,
                temperature=self.temperature,
                main_phases_compositions=self.main_phases_compositions[key],
                reference_phase_compositions=self.reference_phase_compositions[
                    key
                ],
                main_phases_molar_fractions=self.main_phases_molar_fractions[
                    key
                ],
                pressures=self.pressures[key],
                alphas=self.alphas[key],
                iterations=self.iterations[key],
                specified_variable=self.specified_variable[key],
                critical_pressures=self.critical_pressures,
                critical_alphas=self.critical_alphas,
            )
        elif key == "alpha" or key == "a":
            return self.alphas
        elif key == "P":
            return self.pressures
        elif key == "z":
            return self.global_composition
        elif key == "x":
            return self.main_phases_compositions
        elif key == "w":
            return self.reference_phase_compositions

    def __repr__(self):
        display(self.df)
        return ""

    def __len__(self):
        return len(self.alphas)


class TXEnvelope:
    """TXEnvelope.

    This class represents a temperature-composition envelope.

    Attributes
    ----------
    pressure : float
        The pressure at which the envelope is defined. [bar]
    global_composition_0 : np.ndarray
        The global composition at the point where :math:`\alpha = 0`.
    global_composition_i : np.ndarray
        The global composition at the point where :math:`\alpha = 1`.
    main_phases_compositions : np.ndarray
        The compositions of the main phases along the envelope.
        Shape is (n_points, n_phases, n_components).
    reference_phase_compositions : np.ndarray
        The compositions of the reference phase along the envelope.
        Shape is (n_points, n_components).
    main_phases_molar_fractions : np.ndarray
        The molar fractions of the main phases along the envelope.
        Shape is (n_points, n_phases).
    temperatures : np.ndarray
        The temperatures along the envelope. [K]
    alphas : np.ndarray
        The molar fraction of the `global_composition_i`, :math:`\alpha`.
        Shape is (n_points,).
    iterations : np.ndarray
        The number of iterations taken to compute the envelope at each point.
        Shape is (n_points,).
    specified_variable : np.ndarray
        The specified variable used to compute the envelope at each point.
        Shape is (n_points,).
    temperature : float
        The temperatures along envelope. [K]
    critical_temperatures : np.ndarray
        The critical temperatures of the envelope. Shape is (n_critical_points,).
    critical_alphas : np.ndarray
        The molar fractions of the `global_composition_i` at the critical
        points, :math:`\alpha`.  Shape is (n_critical_points,).
    df : pd.DataFrame
        A DataFrame containing the data of the envelope. The columns are:
        - 'alpha': Molar fraction of the `global_composition_i`.
        - 'T': Temperatures along the envelope.
        - 'x_i^j': Compositions of the main phases, where `i`
        is the component index and `j` is the phase index.
        - 'w_i': Compositions of the reference phase, where `i` is the
        component index.
        - 'beta^j': Molar fractions of the main phases, where `j`
        is the phase index.
    """

    def __init__(
        self,
        pressure,
        global_composition_0,
        global_composition_i,
        main_phases_compositions,
        reference_phase_compositions,
        main_phases_molar_fractions,
        temperatures,
        alphas,
        iterations,
        specified_variable,
        critical_temperatures,
        critical_alphas,
        main_phases_kinds,
        reference_phase_kinds,
    ):

        msk = ~np.isnan(temperatures)
        msk_cp = ~np.isnan(critical_temperatures)

        self.temperature = pressure
        self.number_of_components = len(global_composition_0)
        self.number_of_phases = main_phases_compositions.shape[1]
        self.global_composition_0 = global_composition_0
        self.global_composition_i = global_composition_i
        self.main_phases_compositions = main_phases_compositions[msk, :, :]
        self.reference_phase_compositions = reference_phase_compositions[
            msk, :
        ]
        self.main_phases_molar_fractions = main_phases_molar_fractions[msk]
        self.temperatures = temperatures[msk]
        self.alphas = alphas[msk]
        self.iterations = iterations[msk]
        self.specified_variable = specified_variable[msk]
        self.critical_temperatures = critical_temperatures[msk_cp]
        self.critical_alphas = critical_alphas[msk_cp]
        self.main_phases_kinds = main_phases_kinds[msk]
        self.reference_phase_kinds = reference_phase_kinds[msk]

        df = pd.DataFrame()

        df["alpha"] = self.alphas
        df["T"] = self.temperatures

        for i in range(self.number_of_components):
            for j in range(self.number_of_phases):
                df[f"x_{i+1}^{j+1}"] = self.main_phases_compositions[:, j, i]
            df[f"w_{i+1}"] = self.reference_phase_compositions[:, i]

        for i in range(self.number_of_phases):
            df[f"beta^{i+1}"] = self.main_phases_molar_fractions[:, i]

        self.df = df

    def plot(self, **plot_kwargs):
        """Plot the envelope.

        Plot the envelope in a matplotlib axis. Using the natural variables of
        the envelope (:math:`alpha` and temperature)  If the `ax` keyword
        argument is provided, it will be used as the axis to plot on.
        Otherwise, the current axis will be used.

        Parameters
        ----------
        plot_kwargs : dict
            Keyword arguments to pass to the `plot` method of the axis.
            This can include line style, color, etc.
        """
        if "ax" in plot_kwargs:
            ax = plot_kwargs["ax"]
            del plot_kwargs["ax"]
        else:
            ax = plt.gca()
        ax.plot(self.alphas, self.temperatures, **plot_kwargs)
        ax.set_xlabel(r"$\alpha$")
        ax.set_ylabel("Temperature [K]")
        for cp in self.cp:
            ax.scatter(self.alphas[cp], self.temperatures[cp], color="black")

    def __getitem__(self, key):
        if "key" in self.__dict__:
            return self.__dict__["key"]
        elif isinstance(key, np.ndarray) or isinstance(key, list):
            return TXEnvelope(
                global_composition_0=self.global_composition_0,
                global_composition_i=self.global_composition_i,
                pressure=self.pressure,
                main_phases_compositions=self.main_phases_compositions[key],
                reference_phase_compositions=self.reference_phase_compositions[
                    key
                ],
                main_phases_molar_fractions=self.main_phases_molar_fractions[
                    key
                ],
                temperatures=self.temperatures[key],
                alphas=self.alphas[key],
                iterations=self.iterations[key],
                specified_variable=self.specified_variable[key],
                critical_temperatures=self.critical_temperatures,
                critical_alphas=self.critical_alphas,
            )
        elif key == "alpha" or key == "a":
            return self.alphas
        elif key == "T":
            return self.temperatures
        elif key == "z":
            return self.global_composition
        elif key == "x":
            return self.main_phases_compositions
        elif key == "w":
            return self.reference_phase_compositions

    def __repr__(self):
        display(self.df)
        return ""

    def __len__(self):
        return len(self.alphas)
