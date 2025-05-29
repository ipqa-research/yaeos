"""UNIFAC Dortmund Module."""

from typing import List

from ugropy import dortmund, writers

from yaeos.core import GeModel
from yaeos.lib import yaeos_c
from yaeos.models.groups import groups_from_dicts


class UNIFACDortmund(GeModel):
    """UNIFAC Dortmund model.

    Please refer to the `yaeos` user documentation for an in-depth look at the
    model's information: https://ipqa-research.github.io/yaeos/page/index.html

    Parameters
    ----------
    molecules : list of dict
        List of dicts with the groups and their amounts for each molecule.

    Example
    -------
    .. code-block:: python

        from yaeos import UNIFACDortmund

        # Groups for water and ethanol
        water = {"H2O": 1}
        ethanol = {"CH3": 1, "CH2": 1, "OH": 1}

        groups = [water, ethanol]

        model = UNIFACDortmund(groups)

        model.ln_gamma([0.5, 0.5], 298.15)
    """

    def __init__(self, molecules: List[dict]) -> None:
        groups = [writers.to_thermo(m, dortmund) for m in molecules]

        (number_of_groups, groups_ids, groups_ammounts) = groups_from_dicts(
            groups
        )
        self.id = yaeos_c.unifac_dortmund(
            ngs=number_of_groups, g_ids=groups_ids, g_v=groups_ammounts
        )
