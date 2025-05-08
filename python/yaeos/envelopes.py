"""Envelopes.

This module contains the classes that wrapp the data structures used to
represent different kinds of phase envelopes.
"""

import numpy as np

import pandas as pd


class PTEnvelope:
    """PTEnvelope.

    This class represents a pressure-temperature envelope.
    """

    def __init__(
        self,
        global_composition,
        main_phases_compositions,
        reference_phase_compositions,
        main_phases_molar_fractions,
        pressures,
        temperatures,
        iterations,
        specified_variable,
    ):

        msk = ~np.isnan(pressures)
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

        idx = []
        for phase in range(self.number_of_phases):
            cp_idx = []
            lnK = np.log(
                self.reference_phase_compositions
                / self.main_phases_compositions[:, phase]
            )

            for i in range(1, len(lnK)):
                if all(lnK[i, :] * lnK[i - 1, :] < 0):
                    cp_idx.append(i)

            idx.append(cp_idx)

        self.cp = idx

        # lnKs1 = np.log(c["w"][1:]/c["x_l"][1:, 1])
        # lnKs0 = np.log(c["w"][1:]/c["x_l"][1:, 0])
        # for i, lnK in enumerate(np.log(c["w"][1:]/c["x_l"][1:, 1])):
        #     lnK2 = np.log(c["w"]/c["x_l"][i, 1])
        #     crit = (lnK * lnK2 < 0).all()
        #     if crit:
        #         print("asd")

    def __getitem__(self, key):
        if "key" in self.__dict__:
            return self.__dict__["key"]
        elif isinstance(key, np.ndarray):
            return PTEnvelope(
                global_composition=self.global_composition,
                main_phases_compositions=self.main_phases_compositions[key],
                reference_phase_compositions=self.reference_phase_compositions[
                    key
                ],
                main_phases_molar_fractions=self.main_phases_molar_fractions[
                    key
                ],
                pressures=self.pressures[key],
                temperatures=self.temperatures[key],
                iterations=self.iterations[key],
                specified_variable=self.specified_variable[key],
            )
        elif key == "T":
            return self.temperatures
        elif key == "P":
            return self.pressures
        elif key == "z":
            return self.global_composition
        elif key == "x":
            return self.main_phases_compositions
        elif key == "w":
            return self.reference_phase_compositions

    def __repr__(self):
        return self.df.__repr__()

    def __mul__(self, other):
        return PTEnvelope(
            global_composition=self.global_composition,
            main_phases_compositions=self.main_phases_compositions * other,
            reference_phase_compositions=self.reference_phase_compositions
            * other,
            main_phases_molar_fractions=self.main_phases_molar_fractions
            * other,
            pressures=self.pressures * other,
            temperatures=self.temperatures * other,
            iterations=self.iterations,
            specified_variable=self.specified_variable,
        )

    def __len__(self):
        return len(self.temperatures)
