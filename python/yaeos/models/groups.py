"""Handling groups for Group Contribution Models."""

from typing import List

import numpy as np


def groups_from_dicts(molecules: List[dict]):
    """Convert list of dicts with groups and their to C-API format.

    Convert list of dicts with groups and their amounts to the format required
    by the C-API.

    Parameters
    ----------
    molecules : list of dict

    Returns
    -------
    number_of_groups : list of int
        Number of groups for each molecule.
    groups_ids : np.ndarray
        Groups ids for each molecule.
    groups_ammounts : np.ndarray
        Groups amounts for each molecule.
    """
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

    return number_of_groups, groups_ids, groups_ammounts
