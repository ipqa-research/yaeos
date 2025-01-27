import numpy as np

from yaeos import PengRobinson76, QMR


def test_compresibility_factor():
    r = 0.08314462618

    tc = [425.12, 304.21]
    pc = [37.96, 73.83]
    w = [0.200164, 0.223621]
    kij = [[0.0, 0.130], [0.130, 0.0]]

    mixing = QMR(kij, np.zeros((2, 2)))

    eos = PengRobinson76(tc, pc, w, mixing)

    # =========================================================================
    # Composition 0.9
    # =========================================================================
    # 100 F
    t = 310.928
    n = [0.9, 0.1]
    p1 = [41.3685, 68.9476, 137.8951, 206.8427, 275.7903]
    z1 = [0.151, 0.248, 0.482, 0.707, 0.926]

    for p, z in zip(p1, z1):
        v = eos.volume(n, p, t)
        z_calc = p * v / (r * t)

        assert np.allclose(z, z_calc, atol=1e-3)

    # 280 F
    t = 410.928
    n = [0.9, 0.1]
    p2 = [68.9476, 137.8951, 206.8427, 275.7903]
    z2 = [0.289, 0.482, 0.665, 0.840]

    for p, z in zip(p2, z2):
        v = eos.volume(n, p, t)
        z_calc = p * v / (r * t)

        assert np.allclose(z, z_calc, atol=1e-3)

    # 460 F
    t = 510.928
    n = [0.9, 0.1]
    p3 = [41.3685, 68.9476, 137.8951, 206.8427, 275.7903]
    z3 = [0.804, 0.696, 0.643, 0.744, 0.869]

    for p, z in zip(p3, z3):
        v = eos.volume(n, p, t)
        z_calc = p * v / (r * t)

        assert np.allclose(z, z_calc, atol=1e-3)

    # =========================================================================
    # Composition 0.5
    # =========================================================================
    # 100 F
    t = 310.928
    n = [0.5, 0.5]
    p4 = [68.9476, 137.8951, 206.8427, 275.7903]
    z4 = [0.215, 0.404, 0.580, 0.750]

    for p, z in zip(p4, z4):
        v = eos.volume(n, p, t)
        z_calc = p * v / (r * t)

        assert np.allclose(z, z_calc, atol=1e-3)

    # 4100 F
    t = 410.928
    n = [0.5, 0.5]
    p5 = [41.3685, 68.9476, 137.8951, 206.8427, 275.7903]
    z5 = [0.782, 0.638, 0.545, 0.645, 0.765]

    for p, z in zip(p5, z5):
        v = eos.volume(n, p, t)
        z_calc = p * v / (r * t)

        assert np.allclose(z, z_calc, atol=1e-3)

    # 510 F
    t = 510.928
    n = [0.5, 0.5]
    p6 = [41.3685, 68.9476, 137.8951, 206.8427, 275.7903]
    z6 = [0.920, 0.870, 0.796, 0.806, 0.877]

    for p, z in zip(p6, z6):
        v = eos.volume(n, p, t)
        z_calc = p * v / (r * t)

        assert np.allclose(z, z_calc, atol=2e-2)


def test_fugacity():
    r = 0.08314462618

    t = 100  # K
    p = 4.119  # bar
    z_v = [0.958, 1.0 - 0.958]
    z_l = [0.5, 0.5]

    model = PengRobinson76([126.1, 190.6], [33.94, 46.04], [0.040, 0.011])

    # Elliot Z value of vapor
    v_vapor = model.volume(z_v, p, t, root="vapor")

    z_comp_v = p * v_vapor / (r * t)

    assert np.allclose(0.9059, z_comp_v, atol=1e-4)

    # Elliot vapor fugacities
    phi_v = np.exp(model.lnphi_pt(z_v, p, t, root="vapor"))

    assert np.allclose(phi_v[0], 0.9162, atol=1e-4)
    assert np.allclose(phi_v[1], 0.8473, atol=1e-4)

    # Elliot liquid fugacities
    phi_l = np.exp(model.lnphi_pt(z_l, p, t, root="liquid"))

    assert np.allclose(phi_l[0], 1.791, atol=1e-3)
    assert np.allclose(phi_l[1], 0.0937, atol=1e-4)


def test_envelope():
    # TODO: test all the envelope
    model = PengRobinson76(
        [305.32, 540.2], [48.72, 27.4], [0.099493, 0.349469]
    )

    z = [0.1, 0.9]
    envelope = model.phase_envelope_pt(z, t0=230.0, p0=5.0)

    assert np.allclose(envelope["Tc"][0], 532.9, atol=1.0)
    assert np.allclose(envelope["Pc"][0], 32.95, atol=0.2)


def test_flash_pt():
    x_obj = [0.32424471950363210, 0.67575531029866709]
    y_obj = [0.91683466155334536, 8.3165368249135715e-002]
    vx = 8.4918883298198036e-002
    vy = 0.32922132295944545

    model = PengRobinson76(
        [190.564, 425.12], [45.99, 37.96], [0.0115478, 0.200164]
    )

    n = [0.4, 0.6]

    flash_res = model.flash_pt(n, 60.0, 294.0)

    assert np.allclose(flash_res["x"], x_obj, atol=1e-4)
    assert np.allclose(flash_res["y"], y_obj, atol=1e-4)
    assert np.allclose(flash_res["Vx"], vx, atol=1e-4)
    assert np.allclose(flash_res["Vy"], vy, atol=1e-4)
