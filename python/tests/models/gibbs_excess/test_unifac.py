import numpy as np

from yaeos import UNIFACVLE


def test_against_caleb_thermo():
    t = 150
    n = [20.0, 70.0, 10.0]

    groups = [{1: 2}, {1: 1, 2: 1, 14: 1}, {28: 1}]

    model = UNIFACVLE(groups)

    ln_gammas_expected = [0.84433781, -0.19063836, -2.93925506]
    ln_gammas = model.ln_gamma(n, t)

    assert np.allclose(ln_gammas, ln_gammas_expected, atol=1e-8)
