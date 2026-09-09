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


def test_temperature_dependence():
    bij20 = np.array(
        [
            [0.0, 43.552026684128634, 97.82405844415928],
            [-48.846029213395745, 0.0, 180.94738666155268],
            [-384.11635874542793, -208.59497463051014, 0.0],
        ]
    )

    bij30 = np.array(
        [
            [0.0, 21.657910812761273, 99.58458376639004],
            [-37.74765519959855, 0.0, 186.155606345583],
            [-379.32149369269956, -233.34490510114676, 0.0],
        ]
    )

    aij = np.array(
        [
            [0.0, 2.4094201446651944, 0.4861465075816882],
            [-1.4009801734684584, 0.0, 0.710500847827416],
            [-3.0410597123328746, 0.9936949460465081, 0.0],
        ]
    )

    dij = np.array(
        [
            [0.0, -0.007712278604397001, -0.0005200301448241018],
            [0.004210661705262751, 0.0, -0.0003180930394505843],
            [0.005903984936450968, -0.005817018267835272, 0.0],
        ]
    )

    rs = np.array([1.972, 10.496, 31.764])
    qs = np.array([2.105, 12.746, 39.178])

    uni20 = UNIQUAC(qs, rs, bij=bij20)
    uni30 = UNIQUAC(qs, rs, bij=bij30)
    uni_linear = UNIQUAC(qs, rs, aij=aij, dij=dij)

    n = np.array([5.0, 15.0, 65.0])

    ln_gamma20 = uni20.ln_gamma(n, 293.15)
    ln_gamma30 = uni30.ln_gamma(n, 303.15)

    ln_gamma_linear20 = uni_linear.ln_gamma(n, 293.15)
    ln_gamma_linear30 = uni_linear.ln_gamma(n, 303.15)

    assert np.allclose(ln_gamma20, ln_gamma_linear20)
    assert np.allclose(ln_gamma30, ln_gamma_linear30)
