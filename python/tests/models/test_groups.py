from yaeos.models.groups import groups_from_dicts


def test_groups_from_dicts():
    molecules = [
        {"1": 1, "2": 2},
        {"1": 1, "2": 2, "3": 3},
        {"1": 1},
    ]
    number_of_groups, groups_ids, groups_ammounts = groups_from_dicts(
        molecules
    )
    assert number_of_groups == [2, 3, 1]
    assert groups_ids.tolist() == [[1, 2, 0], [1, 2, 3], [1, 0, 0]]
    assert groups_ammounts.tolist() == [[1, 2, 0], [1, 2, 3], [1, 0, 0]]
