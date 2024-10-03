import numpy as np

from yaeos import UNIQUAC


def test_ln_gammas():
    a = np.array(
        [
            [0.0, -75.46, -60.15],
            [120.20, 0.0, 44.22],
            [120.20, 33.21, 0.0],
        ]
    )

    b = np.array(
        [
            [0.0, -0.10062, 0.2566],
            [0.44835, 0.0, -0.01325],
            [0.44835, 0.124, 0.0],
        ]
    )

    c = np.array(
        [
            [0.0, -0.0008052, 0.00021],
            [0.0004704, 0.0, -0.00033],
            [0.0004704, -0.000247, 0.0],
        ]
    )

    d = np.array(
        [
            [0.0, -0.001, 0.0002],
            [-0.001, 0.0, 0.0002],
            [-0.001, 0.0002, 0.0],
        ]
    )

    e = np.array(
        [
            [0.0, -0.00001, 0.00001],
            [-0.00001, 0.0, 0.00001],
            [-0.00001, 0.00001, 0.0],
        ]
    )

    rs = np.array([0.92, 2.1055, 1.5])
    qs = np.array([1.4, 1.972, 1.4])

    model = UNIQUAC(qs, rs, a, b, c, d, e)

    t = 298.15
    n = np.array([20.0, 70.0, 10.0])

    ln_gammas = model.ln_gamma(n, t)

    expect = np.array(
        [-164.62277497059728, -60.906444787104235, -75.52457152449654]
    )

    assert np.allclose(ln_gammas, expect)


def test_give_only_bij():
    b = np.array(
        [
            [0.0, -526.02, -309.64],
            [318.06, 0.0, 91.532],
            [-1325.1, -302.57, 0.0],
        ]
    )

    rs = np.array([0.92, 2.1055, 3.1878])
    qs = np.array([1.4, 1.972, 2.4])

    t = 298.15

    model = UNIQUAC(qs, rs, bij=b)

    n = np.array([2.0, 2.0, 8.0])

    gammas = np.exp(model.ln_gamma(n, t))

    excepted = np.array([8.856, 0.860, 1.425])

    assert np.allclose(gammas, excepted, atol=1e-3)
