import numpy as np
import numpy.testing as npt

from yaeos import SoaveRedlichKwong, PengRobinson76, QMR


def test_dummy():
    mr = QMR(np.zeros((2, 2)), np.zeros((2, 2)))

    model = PengRobinson76(
        np.array([300, 350]), np.array([30, 40]), np.array([0.152, 0.325]), mr
    )

    fug = model.lnphi_vt(np.array([5, 6]), 2.0, 303.15)["ln_phi"]

    assert np.allclose(fug, np.array([0.65960035, 0.30581106]))


def test_flash():
    nc = 2
    n = [0.4, 0.6]
    Tc = [190.564, 425.12]
    Pc = [45.99, 37.96]
    w = [0.0115478, 0.200164]

    lij = kij = np.zeros((nc, nc))
    mixrule = QMR(kij, lij)
    model = PengRobinson76(Tc, Pc, w, mixrule)

    P, T = 60.0, 294.0
    flash = model.flash_pt(n, pressure=P, temperature=T)

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
