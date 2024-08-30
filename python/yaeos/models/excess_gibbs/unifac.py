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
        groups_ids = np.zeros((nc, max_ng), dtype=np.int32, order="F")
        groups_ammounts = np.zeros((nc, max_ng), dtype=np.int32, order="F")

        ngs = []
        for i, grp in enumerate(molecules):
            ids = []
            vs = []
            for group_id, group_count in grp.items():
                ids.append(int(group_id))
                vs.append(int(group_count))

            ngs.append(len(ids))
            groups_ids[i, : len(ids)] = ids
            groups_ammounts[i, : len(ids)] = vs

        self.id = yaeos_c.unifac_vle(
            ngs=ngs, g_ids=groups_ids, g_v=groups_ammounts
        )
