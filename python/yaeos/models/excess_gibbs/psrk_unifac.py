"""PSRK-UNIFAC Module."""

from typing import List

from yaeos.core import GeModel
from yaeos.lib import yaeos_c
from yaeos.models.groups import groups_from_dicts


class UNIFACPSRK(GeModel):
    """PSRK-UNIFAC model.

    Please refer to the `yaeos` user documentation for an in-depth look at the
    model's information: https://ipqa-research.github.io/yaeos/page/index.html

    Parameters
    ----------
    molecules : list of dict
        List of dicts with the groups and their amounts for each molecule.

    Example
    -------
    .. code-block:: python

        from yaeos import UNIFACPSRK

        # Groups for water and ethanol
        water = {16: 1}
        ethanol = {1: 1, 2: 1, 14: 1}

        groups = [water, ethanol]

        model = UNIFACPSRK(groups)

        model.ln_gamma([0.5, 0.5], 298.15)
    """

    def __init__(self, molecules: List[dict]) -> None:

        (number_of_groups, groups_ids, groups_ammounts) = groups_from_dicts(
            molecules
        )
        self.id = yaeos_c.unifac_psrk(
            ngs=number_of_groups, g_ids=groups_ids, g_v=groups_ammounts
        )
