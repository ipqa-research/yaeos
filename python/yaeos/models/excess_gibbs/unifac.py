"""UNIFAC Module.
"""

import numpy as np

from yaeos.core import GeModel
from yaeos.lib import yaeos_c


class UNIFACVLE(GeModel):
    """UNIFAC VLE model.

    Parameters
    ----------
    molecules : list of dict
        List of dicts with the groups and their amounts for each molecule.
    """

    def __init__(self, molecules) -> None:

        nc = len(molecules)
        max_ng = max([len(groups) for groups in molecules])

        # The C-API expects two 2D arrays with the groups and their amounts
        # for each molecule. Each row in the arrays corresponds to a molecule
        # and the columns to the groups. The arrays are padded with zeros.
        groups_ids = np.zeros((nc, max_ng), dtype=np.int32, order="F")
        groups_ammounts = np.zeros((nc, max_ng), dtype=np.int32, order="F")

        # Construction of the arrays from the input dictionary
        number_of_groups = []
        for i, molecule_groups in enumerate(molecules):
            ids = []
            ammount = []
            for group_id, group_count in molecule_groups.items():
                ids.append(int(group_id))
                ammount.append(int(group_count))

            number_of_groups.append(len(ids))
            groups_ids[i, : len(ids)] = ids
            groups_ammounts[i, : len(ids)] = ammount

        self.id = yaeos_c.unifac_vle(
            ngs=number_of_groups, g_ids=groups_ids, g_v=groups_ammounts
        )
