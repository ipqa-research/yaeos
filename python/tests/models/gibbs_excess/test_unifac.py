import numpy as np

from yaeos import UNIFACVLE


def test_against_caleb_thermo():
    t = 150
    n = [20.0, 70.0, 10.0]

    nt = np.sum(n)

    groups = [{1: 2}, {1: 1, 2: 1, 14: 1}, {28: 1}]

    model = UNIFACVLE(groups)

    # =========================================================================
    # ln_gammas
    # =========================================================================
    ln_gammas_expected = [
        0.8443378123549767,
        -0.19063835833404139,
        -2.9392550622509956,
    ]
    ln_gammas = model.ln_gamma(n, t)

    assert np.allclose(ln_gammas, ln_gammas_expected, atol=1e-8)

    # =========================================================================
    # Ge
    # =========================================================================
    ge, deriv = model.excess_gibbs(
        n, t, dt=True, dt2=True, dn=True, dtn=True, dn2=True
    )

    gd1, d1 = model.excess_gibbs(n, t, dt=True)
    gd2, d2 = model.excess_gibbs(n, t, dt2=True)
    gd3, d3 = model.excess_gibbs(n, t, dn=True)
    gd4, d4 = model.excess_gibbs(n, t, dtn=True)
    gd5, d5 = model.excess_gibbs(n, t, dn2=True)

    assert np.isclose(ge, gd1)
    assert np.isclose(ge, gd2)
    assert np.isclose(ge, gd3)
    assert np.isclose(ge, gd4)
    assert np.isclose(ge, gd5)

    assert np.isclose(deriv["dt"], d1["dt"])
    assert np.isclose(deriv["dt2"], d2["dt2"])
    assert np.allclose(deriv["dn"], d3["dn"])
    assert np.allclose(deriv["dtn"], d4["dtn"])
    assert np.allclose(deriv["dn2"], d5["dn2"])

    assert np.isclose(-3.223992676822129, ge / nt)
    assert np.isclose(0.03268447167877294, deriv["dt"] / nt)
    assert np.isclose(-0.0003594405355829625, deriv["dt2"] / nt)

    r = 0.08314462618
    assert np.allclose(ln_gammas_expected, deriv["dn"] / r / t)

    assert np.allclose(
        [0.06015388937032817, 0.02239721670764982, 0.0497564210935243],
        deriv["dtn"],
    )

    dgedn2_exp = np.array(
        [
            [-9.38494059, 1.67630836, 7.03572265],
            [1.67630836, 4.32872363, -33.65368214],
            [7.03572265, -33.65368214, 221.50432968],
        ]
    )

    assert np.allclose(dgedn2_exp, deriv["dn2"] * nt)

    # =========================================================================
    # He
    # =========================================================================
    he, deriv = model.excess_enthalpy(n, t, dt=True, dn=True)

    hd1, d1 = model.excess_enthalpy(n, t, dt=True)
    hd2, d2 = model.excess_enthalpy(n, t, dn=True)

    assert np.isclose(he, hd1)
    assert np.isclose(he, hd2)

    assert np.isclose(deriv["dt"], d1["dt"])
    assert np.allclose(deriv["dn"], d2["dn"])

    assert np.isclose(-8.126663428638069, he / nt)
    assert np.isclose(0.05391608033744438, deriv["dt"] / nt)
    assert np.allclose(
        [1.5072393613288813, -5.737165762079205, -44.120952674484016],
        deriv["dn"],
    )

    # =========================================================================
    # Se
    # =========================================================================
    se, deriv = model.excess_entropy(n, t, dt=True, dn=True)

    sd1, d1 = model.excess_entropy(n, t, dt=True)
    sd2, d2 = model.excess_entropy(n, t, dn=True)

    assert np.isclose(se, sd1)
    assert np.isclose(se, sd2)

    assert np.isclose(deriv["dt"], d1["dt"])
    assert np.allclose(deriv["dn"], d2["dn"])

    assert np.isclose(-0.03268447167877293, se / nt)
    assert np.isclose(0.0003594405355829625, deriv["dt"] / nt)
    assert np.allclose(
        [-0.06015388937032817, -0.02239721670764982, -0.0497564210935243],
        deriv["dn"],
    )

    # =========================================================================
    # Cpe
    # =========================================================================
    cpe = model.excess_cp(n, t)

    assert np.isclose(0.05391608033744438, cpe / nt)
