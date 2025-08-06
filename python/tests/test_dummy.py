import numpy as np
import numpy.testing as npt

from yaeos import PengRobinson76, QMR, SoaveRedlichKwong


def test_dummy():
    mr = QMR(np.zeros((2, 2)), np.zeros((2, 2)))

    model = PengRobinson76(
        np.array([300, 350]), np.array([30, 40]), np.array([0.152, 0.325]), mr
    )

    fug = model.lnphi_vt(np.array([5, 6]), 2.0, 303.15)

    assert np.allclose(fug, np.array([0.65960035, 0.30581106]))


def test_flash():
    nc = 2
    n = [0.4, 0.6]
    tc = [190.564, 425.12]
    pc = [45.99, 37.96]
    w = [0.0115478, 0.200164]

    lij = kij = np.zeros((nc, nc))
    mixrule = QMR(kij, lij)
    model = PengRobinson76(tc, pc, w, mixrule)

    p, t = 60.0, 294.0
    flash = model.flash_pt(n, pressure=p, temperature=t)

    npt.assert_allclose(flash["x"], [0.32424472, 0.67575528], rtol=1e-5)


def test_saturation():
    nc = 2
    n = [0.4, 0.6]
    tc = [190.564, 425.12]
    pc = [45.99, 37.96]
    w = [0.0115478, 0.200164]

    lij = kij = np.zeros((nc, nc))
    mixrule = QMR(kij, lij)
    model = SoaveRedlichKwong(tc, pc, w, mixrule)

    bub = model.saturation_pressure(n, temperature=303.15, kind="bubble")
    npt.assert_allclose(bub["y"], [0.899642, 0.100358], rtol=1e-5)


def test_envelope_pt():
    import yaeos
    import numpy as np

    tc = [304.21, 373.53, 190.564]
    pc = [73.83000000000001, 89.62910000000001, 45.99]
    w = [0.223621, 0.0941677, 0.0115478]

    kij = np.zeros((3, 3))

    kij[0, 1] = kij[1, 0] = 0.0974
    kij[0, 2] = kij[2, 0] = 0.110
    kij[1, 2] = kij[2, 1] = 0.069

    mr = yaeos.QMR(kij, 0 * kij)

    model = yaeos.PengRobinson76(tc, pc, w, mr)
    z = np.array([0.0987, 0.4023, 0.4990])
    z = z / sum(z)
    dew = model.phase_envelope_pt(z, kind="dew", p0=0.01, max_points=1500)

    assert dew["P"][-1] > 1000
    npt.assert_allclose(dew["Tc"], [275, 218], atol=1)
    npt.assert_allclose(dew["Pc"], [122, 95], atol=1)
