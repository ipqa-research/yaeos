"""Gerg2008 Equation of State"""

from yaeos.core import ArModel
from yaeos.lib import yaeos_c


class GERG2008(ArModel):
    possible_components = {
        "methane": 1,
        "c1": 1,
        "ch4": 1,
        "nitrogen": 2,
        "n2": 2,
        "carbon dioxide": 3,
        "co2": 3,
        "ethane": 4,
        "c2": 4,
        "propane": 5,
        "c3": 5,
        "n-butane": 6,
        "nbutane": 6,
        "nc4": 6,
        "isobutane": 7,
        "ic4": 7,
        "i-butane": 7,
        "n-pentane": 8,
        "npentane": 8,
        "nc5": 8,
        "isopentane": 9,
        "ic5": 9,
        "i-pentane": 9,
        "n-hexane": 10,
        "nhexane": 10,
        "nc6": 10,
        "n-heptane": 11,
        "nheptane": 11,
        "nc7": 11,
        "n-octane": 12,
        "noctane": 12,
        "nc8": 12,
        "nonane": 13,
        "c9": 13,
        "nc9": 13,
        "decane": 14,
        "c10": 14,
        "nc10": 14,
        "hydrogen": 15,
        "h2": 15,
        "oxygen": 16,
        "o2": 16,
        "carbon monoxide": 17,
        "co": 17,
        "water": 18,
        "h2o": 18,
        "hydrogen sulfide": 19,
        "h2s": 19,
        "helium": 20,
        "he": 20,
        "argon": 21,
        "ar": 21,
    }

    def __init__(self, names):
        """GERG2008 Equation of State.

        Parameters
        ----------
        names: list, str
            List of names for the components to use in the model.
        """

        ids = [self.possible_components[name] for name in names]
        self.id = yaeos_c.multifluid_gerg2008(ids)
