import pytest

from yaeos import QMR


def test_making_it_explode():
    mr = QMR([[0.0, 0.1], [0.1, 0.0]], [[0.0, 0.1], [0.1, 0.0]])

    with pytest.raises(NotImplementedError):
        super(QMR, mr).set_mixrule(0)
