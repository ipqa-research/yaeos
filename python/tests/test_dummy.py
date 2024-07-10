import numpy as np

from yaeos import PengRobinson76, QMR


def test_dummy():
    mr = QMR(np.zeros((2, 2)), np.zeros((2, 2)))

    model = PengRobinson76(
        np.array([300, 350]), np.array([30, 40]), np.array([0.152, 0.325]), mr
    )

    fug = model.fugacity(np.array([5, 6]), 2.0, 303.15)["ln_phi"]

    assert np.allclose(fug, np.array([2.84731863, 2.49352934]))
