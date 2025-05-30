import numpy as np

from yaeos import PSRK


def test_psrk_paper2005_1():
    tc = [456.8, 647.3]
    pc = [53.9049, 220.48321]
    w = [0.407, 0.344]

    c1 = [1.62716, 1.0783]
    c2 = [-5.00731, -0.58321]
    c3 = [10.34189, 0.54619]

    grps = [{145: 1}, {16: 1}]

    model = PSRK(tc, pc, w, grps, c1, c2, c3)

    x1_exp = [
        0.00754,
        0.04774,
        0.07538,
        0.12312,
        0.17839,
        0.26382,
        0.33668,
        0.4196,
        0.5402,
        0.63317,
        0.74623,
        0.86935,
        0.94221,
        0.97739,
    ]

    pbub_exp = [
        0.06384039999999999,
        0.2374065,
        0.32718200000000003,
        0.4349127,
        0.5067332,
        0.5645884999999999,
        0.5865337,
        0.5985037,
        0.6104738,
        0.6244389,
        0.6503741,
        0.6962594,
        0.7341646000000001,
        0.7541146999999999,
    ]

    y1_exp = [
        0.01759,
        0.19849,
        0.3593,
        0.5,
        0.57789,
        0.68844,
        0.76382,
        0.83417,
        0.88442,
        0.91457,
        0.94472,
        0.95729,
        0.96482,
        0.97236,
        0.9799,
        0.98995,
    ]

    pdew_exp = [
        0.021945100000000002,
        0.035910199999999996,
        0.0478803,
        0.0578554,
        0.06384039999999999,
        0.07182039999999999,
        0.09177060000000001,
        0.12169579999999999,
        0.1715711,
        0.2354115,
        0.3431421,
        0.42892769999999997,
        0.5087282,
        0.6024938,
        0.6663342,
        0.7261845,
    ]

    t = 291.15

    pbub = []
    for i, x in enumerate(x1_exp):
        sol = model.saturation_pressure([x, 1 - x], t, p0=pbub_exp[i])
        pbub.append(sol["P"])

    assert np.allclose(pbub, pbub_exp, atol=4e-3)

    pdew = []
    for i, y in enumerate(y1_exp):
        sol = model.saturation_pressure(
            [y, 1 - y], t, kind="dew", p0=pdew_exp[i]
        )
        pdew.append(sol["P"])

    assert np.allclose(pdew, pdew_exp, atol=2e-2)


def test_psrk_paper2005_2():
    tc = [324.6, 417.0]
    pc = [83.0865, 77.007]
    w = [0.12, 0.073]

    c1 = [0.66635, 0.55192]
    c2 = [0.35497, 0.01934]
    c3 = [-1.3766, 0.59414]

    grps = [{130: 1}, {143: 1}]

    x1_exp = [
        0.01942,
        0.07864,
        0.14284,
        0.20951,
        0.28609,
        0.36516,
        0.44671,
        0.53322,
        0.62219,
        0.73341,
        0.84957,
        0.96077,
    ]

    pbub_exp = [
        4.6423,
        6.8492,
        8.8347,
        10.746199999999998,
        12.5833,
        14.2728,
        15.814699999999998,
        17.43,
        19.045,
        21.027,
        23.229799999999997,
        25.506600000000002,
    ]

    y1_exp = [
        0.06402,
        0.1432,
        0.23227,
        0.30896,
        0.39059,
        0.45983,
        0.54885,
        0.6131,
        0.67487,
        0.72425,
        0.77361,
        0.818,
        0.86483,
        0.91409,
        0.9535,
        0.9806,
    ]

    pdew_exp = [
        4.049300000000001,
        4.412,
        4.9215,
        5.4318,
        6.089099999999999,
        6.8210999999999995,
        7.9939,
        9.1686,
        10.564499999999999,
        12.0351,
        13.874200000000002,
        15.861099999999999,
        18.2901,
        21.2348,
        23.5907,
        25.136599999999998,
    ]

    model = PSRK(tc, pc, w, grps, c1, c2, c3)

    t = 273.15

    pbub = []
    for i, x in enumerate(x1_exp):
        sol = model.saturation_pressure([x, 1 - x], t, p0=pbub_exp[i])
        pbub.append(sol["P"])

    assert np.allclose(pbub, pbub_exp, rtol=1.1e-2)

    pdew = []
    for i, y in enumerate(y1_exp):
        sol = model.saturation_pressure(
            [y, 1 - y], t, kind="dew", p0=pdew_exp[i]
        )
        pdew.append(sol["P"])

    assert np.allclose(pdew, pdew_exp, rtol=7e-3)


def test_psrk_paper2005_3():
    tc = [324.6, 417.0]
    pc = [83.0865, 77.007]
    w = [0.12, 0.073]

    c1 = [0.66635, 0.55192]
    c2 = [0.35497, 0.01934]
    c3 = [-1.3766, 0.59414]

    grps = [{130: 1}, {143: 1}]

    x1_exp = [
        0.02719,
        0.10883,
        0.17316,
        0.24987,
        0.31174,
        0.36618,
        0.42805,
        0.50477,
        0.58149,
        0.67553,
        0.81659,
        0.89826,
        0.96508,
    ]

    pbub_exp = [
        0.514,
        0.9501999999999999,
        1.2403,
        1.5295,
        1.7460999999999998,
        1.8895,
        2.0322999999999998,
        2.1741,
        2.3896,
        2.62,
        2.8882000000000003,
        3.1033,
        3.2458,
    ]

    y1_exp = [
        0.02968,
        0.10641,
        0.18561,
        0.26977,
        0.35888,
        0.41581,
        0.52719,
        0.61876,
        0.69548,
        0.79446,
        0.86621,
        0.92062,
        0.97006,
    ]

    pdew_exp = [
        0.2927,
        0.287,
        0.3549,
        0.3487,
        0.4158,
        0.4116,
        0.4771,
        0.6178,
        0.7595000000000001,
        1.1208,
        1.5578,
        2.1434,
        2.8769,
    ]

    model = PSRK(tc, pc, w, grps, c1, c2, c3)

    t = 213.15

    pbub = []
    for i, x in enumerate(x1_exp):
        sol = model.saturation_pressure([x, 1 - x], t, p0=pbub_exp[i])
        pbub.append(sol["P"])

    assert np.allclose(pbub, pbub_exp, rtol=1e-1)

    pdew = []
    for i, y in enumerate(y1_exp):
        sol = model.saturation_pressure(
            [y, 1 - y], t, kind="dew", p0=pdew_exp[i]
        )
        pdew.append(sol["P"])

    assert np.allclose(pdew, pdew_exp, rtol=0.2)


def test_psrk_paper2005_4():
    tc = [440.0, 523.2]
    pc = [91.1925, 38.300850000000004]
    w = [0.318, 0.363]

    c1 = [0.96273, 1.0408]
    c2 = [0.10351, -0.17686]
    c3 = [-1.6102, 0.49506]

    grps = [{149: 1}, {1: 1, 2: 1, 21: 1}]

    x1_exp = [
        0.03202,
        0.12562,
        0.24877,
        0.34975,
        0.46798,
        0.58867,
        0.68473,
        0.80049,
        0.94089,
        0.9803,
    ]

    pbub_exp = [
        0.014705900000000001,
        0.0392157,
        0.0784314,
        0.11764709999999999,
        0.17279409999999998,
        0.2389706,
        0.2965686,
        0.372549,
        0.4607843,
        0.48284309999999997,
    ]

    y1_exp = [
        0.03695,
        0.13547,
        0.26601,
        0.45074,
        0.61576,
        0.74877,
        0.85714,
        0.94335,
        0.97291,
        0.98768,
        0.99261,
        1.0,
    ]

    pdew_exp = [
        0.0073529,
        0.0085784,
        0.009803899999999999,
        0.013480399999999998,
        0.019607799999999998,
        0.0281863,
        0.0465686,
        0.09313729999999999,
        0.1446078,
        0.21446079999999998,
        0.2830882,
        0.47058819999999996,
    ]

    model = PSRK(tc, pc, w, grps, c1, c2, c3)

    t = 251.95

    pbub = []
    for i, x in enumerate(x1_exp):
        sol = model.saturation_pressure([x, 1 - x], t, p0=pbub_exp[i])
        pbub.append(sol["P"])

    assert np.allclose(pbub, pbub_exp, rtol=2e-2)

    pdew = []
    for i, y in enumerate(y1_exp):
        sol = model.saturation_pressure(
            [y, 1 - y], t, kind="dew", p0=pdew_exp[i]
        )
        pdew.append(sol["P"])

    assert np.allclose(pdew, pdew_exp, rtol=1)
