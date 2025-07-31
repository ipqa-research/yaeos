from typing import List

import numpy as np

from scipy.optimize import approx_fprime

import yaeos
from yaeos.core import ArModel


def test_same_as_caleb_vt():
    mods = [
        yaeos.PengRobinson76,
        yaeos.PengRobinson78,
        yaeos.SoaveRedlichKwong,
    ]

    tc = [469.7, 507.6]
    pc = [33.7, 30.250000000000004]
    w = [0.251506, 0.301261]

    n = np.array([4.0, 6.0])
    nt = np.sum(n)

    v = 1.0 * nt * 1000  # Liters (1 m^3/mol * nt [mol] * 1000 L/m^3)
    t = 298.15

    models: List[ArModel] = []

    for mod in mods:
        model = mod(
            critical_temperatures=tc, critical_pressures=pc, acentric_factors=w
        )
        models.append(model)

    # =========================================================================
    # Pressure
    # =========================================================================
    p_exp = [0.024758433028831645, 0.024758433028831645, 0.024759435993050086]
    dpdv_exp = [
        -2.472730282235202e-05,
        -2.472730282235202e-05,
        -2.472930571514075e-05,
    ]
    dpdt_exp = [
        8.321187984959119e-05,
        8.321187984959119e-05,
        8.321898105202635e-05,
    ]
    dpdn_exp = [
        [0.024733557386271844, 0.024723133113072143],
        [0.024733557386271844, 0.024723133113072143],
        [0.024735507908033087, 0.02472517091987919],
    ]

    for idx, model in enumerate(models):
        pi = model.pressure(n, v, t)

        pd1, der1 = model.pressure(n, v, t, dv=True)
        pd2, der2 = model.pressure(n, v, t, dt=True)
        pd3, der3 = model.pressure(n, v, t, dn=True)
        pd4, der4 = model.pressure(n, v, t, dv=True, dt=True, dn=True)

        assert np.isclose(pi, pd1)
        assert np.isclose(pi, pd2)
        assert np.isclose(pi, pd3)
        assert np.isclose(pi, pd4)

        assert np.isclose(der1["dv"], der4["dv"])
        assert np.isclose(der2["dt"], der4["dt"])
        assert np.allclose(der3["dn"], der4["dn"])

        assert np.isclose(p_exp[idx], pi)
        assert np.isclose(dpdv_exp[idx] / nt, der4["dv"])  # Caleb is V permole
        assert np.isclose(dpdt_exp[idx], der4["dt"])
        # The next is over nt cause Caleb is a constant molar volume
        assert np.allclose(dpdn_exp[idx] / nt, der4["dn"])

    # =========================================================================
    # Residual enthalpy VT
    # =========================================================================
    h_exp = [-0.08233139437821592, -0.08233139437821592, -0.08244045745725088]
    dhdv_exp = [
        8.231915480358669e-05,
        8.231915480358669e-05,
        8.243348552089947e-05,
    ]
    dhdt_exp = [
        0.00011201776832480803,
        0.00011201776832480803,
        0.00012596855942014428,
    ]

    epsilon = 1e-6

    for idx, model in enumerate(models):
        # Caleb doesn't have dH_dn at constant V
        dhdn_exp = approx_fprime(n, model.enthalpy_residual_vt, epsilon, v, t)

        hi = model.enthalpy_residual_vt(n, v, t)

        hd1, der1 = model.enthalpy_residual_vt(n, v, t, dv=True)
        hd2, der2 = model.enthalpy_residual_vt(n, v, t, dt=True)
        hd3, der3 = model.enthalpy_residual_vt(n, v, t, dn=True)
        hd4, der4 = model.enthalpy_residual_vt(
            n, v, t, dv=True, dt=True, dn=True
        )

        assert np.isclose(hi, hd1)
        assert np.isclose(hi, hd2)
        assert np.isclose(hi, hd3)
        assert np.isclose(hi, hd4)

        assert np.isclose(der1["dv"], der4["dv"])
        assert np.isclose(der2["dt"], der4["dt"])
        assert np.allclose(der3["dn"], der4["dn"])

        assert np.isclose(h_exp[idx], hi / nt, atol=4e-6)
        assert np.isclose(dhdv_exp[idx], der4["dv"])  # Caleb is permole
        assert np.isclose(dhdt_exp[idx], der4["dt"] / nt, atol=2e-8)
        assert np.allclose(dhdn_exp, der4["dn"])

    # =========================================================================
    # Residual entropy VT
    # =========================================================================
    s_exp = [
        -6.725919527969187e-05,
        -6.725919527969187e-05,
        -7.435800385866321e-05,
    ]

    epsilon = 1e-6

    for idx, model in enumerate(models):

        def s_v(v, n, t):
            return model.entropy_residual_vt(n, v[0], t)

        def s_t(t, n, v):
            return model.entropy_residual_vt(n, v, t[0])

        dsdv_exp = approx_fprime([v], s_v, epsilon, n, t)
        dsdt_exp = approx_fprime([t], s_t, epsilon, n, v)
        dsdn_exp = approx_fprime(n, model.entropy_residual_vt, epsilon, v, t)

        si = model.entropy_residual_vt(n, v, t)

        sd1, der1 = model.entropy_residual_vt(n, v, t, dv=True)
        sd2, der2 = model.entropy_residual_vt(n, v, t, dt=True)
        sd3, der3 = model.entropy_residual_vt(n, v, t, dn=True)
        sd4, der4 = model.entropy_residual_vt(
            n, v, t, dv=True, dt=True, dn=True
        )

        assert np.isclose(si, sd1)
        assert np.isclose(si, sd2)
        assert np.isclose(si, sd3)
        assert np.isclose(si, sd4)

        assert np.isclose(der1["dv"], der4["dv"])
        assert np.isclose(der2["dt"], der4["dt"])
        assert np.allclose(der3["dn"], der4["dn"])

        assert np.isclose(s_exp[idx], si / nt)
        assert np.isclose(dsdv_exp, der4["dv"])
        assert np.isclose(dsdt_exp, der4["dt"])
        assert np.allclose(dsdn_exp, der4["dn"])

    # =========================================================================
    # Residual Gibbs VT
    # =========================================================================
    g_exp = [-0.06227806530557579, -0.06227806530557579, -0.06027061860679045]

    epsilon = 1e-6

    for idx, model in enumerate(models):

        def g_v(v, n, t):
            return model.gibbs_residual_vt(n, v[0], t)

        def g_t(t, n, v):
            return model.gibbs_residual_vt(n, v, t[0])

        dgdv_exp = approx_fprime([v], g_v, 1e-4, n, t)
        dgdt_exp = approx_fprime([t], g_t, epsilon, n, v)
        dgdn_exp = approx_fprime(n, model.gibbs_residual_vt, epsilon, v, t)

        gi = model.gibbs_residual_vt(n, v, t)
        si = model.entropy_residual_vt(n, v, t)
        hi = model.enthalpy_residual_vt(n, v, t)

        # just because
        assert np.isclose(gi, hi - t * si)

        gd1, der1 = model.gibbs_residual_vt(n, v, t, dv=True)
        gd2, der2 = model.gibbs_residual_vt(n, v, t, dt=True)
        gd3, der3 = model.gibbs_residual_vt(n, v, t, dn=True)
        gd4, der4 = model.gibbs_residual_vt(n, v, t, dv=True, dt=True, dn=True)

        assert np.isclose(gi, gd1)
        assert np.isclose(gi, gd2)
        assert np.isclose(gi, gd3)
        assert np.isclose(gi, gd4)

        assert np.isclose(der1["dv"], der4["dv"])
        assert np.isclose(der2["dt"], der4["dt"])
        assert np.allclose(der3["dn"], der4["dn"])

        assert np.isclose(g_exp[idx], gi / nt, atol=2e-6)
        assert np.isclose(dgdv_exp, der4["dv"])
        assert np.isclose(dgdt_exp, der4["dt"])
        assert np.allclose(dgdn_exp, der4["dn"])

    # =========================================================================
    # Residual Cp
    # =========================================================================
    cp_exp = [
        0.0003890367245792703,
        0.0003890367245792703,
        0.0004033734630853303,
    ]

    for idx, model in enumerate(models):
        cp = model.cp_residual_vt(n, v, t)

        assert np.isclose(cp_exp[idx], cp / nt, atol=1e-7)

    # =========================================================================
    # Residual Cv
    # =========================================================================
    cv_exp = [
        4.476410026600885e-05,
        4.476410026600885e-05,
        5.161368892615888e-05,
    ]

    for idx, model in enumerate(models):
        cv = model.cv_residual_vt(n, v, t)

        assert np.isclose(cv_exp[idx], cv / nt)

    # =========================================================================
    # ln_phi
    # =========================================================================
    ln_phi_e = [
        [-0.0010030671414005777, -0.0014236487978196851],
        [-0.0010030671414005777, -0.0014236487978196851],
        [-0.00096472785584406, -0.0013817578146076205],
    ]

    for idx, model in enumerate(models):
        ln_phi = model.lnphi_vt(n, v, t)

        assert np.allclose(ln_phi_e[idx], ln_phi, atol=1e-7)
