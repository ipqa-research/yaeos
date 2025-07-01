from pathlib import Path

import numpy as np

import pytest

from yaeos import NRTL, UNIFACDortmund, UNIFACPSRK, UNIFACVLE, UNIQUAC


data_path = (
    Path(__file__).parent.parent.parent.parent.parent / "ge_test_vals.txt"
)


def test_same_as_fortran():
    # =========================================================================
    # Read data
    # =========================================================================
    if not data_path.exists():
        pytest.skip("Test data not found: GeModels same as fortran")

    with open(data_path, "r") as f:
        data_lines = f.readlines()

    # =========================================================================
    # Setup models
    # =========================================================================
    # NRTL
    a = np.array(
        [
            [0.0, -0.801, -0.351],
            [-0.523, 0.0, 0.214],
            [0.127, 0.211, 0.0],
        ]
    )

    b = np.array(
        [
            [0.0, -586.1, 246.2],
            [301.2, 0.0, -104.2],
            [150.23, -114.78, 0.0],
        ]
    )

    c = np.array(
        [
            [0.0, 0.3, 0.3],
            [0.3, 0.0, 0.3],
            [0.3, 0.3, 0.0],
        ]
    )

    nrtl = NRTL(a, b, c)

    # UNIFAC VLE
    groups = [{1: 2, 2: 4}, {1: 1, 2: 1, 14: 1}, {9: 5, 11: 1}]

    unifac = UNIFACVLE(groups)

    # PSRK Dortmund
    psrk = UNIFACPSRK(groups)

    # UNIFAC Dortmund
    unifac_dortmund = UNIFACDortmund(groups)

    # UNIQUAC
    aij = np.array(
        [
            [0.0, -75.46, -60.15],
            [120.20, 0.0, 44.22],
            [120.20, 33.21, 0.0],
        ]
    )

    bij = np.array(
        [
            [0.0, -0.10062, 0.2566],
            [0.44835, 0.0, -0.01325],
            [0.44835, 0.124, 0.0],
        ]
    )

    cij = np.array(
        [
            [0.0, -0.0008052, 0.00021],
            [0.0004704, 0.0, -0.00033],
            [0.0004704, -0.000247, 0.0],
        ]
    )

    dij = np.array(
        [
            [0.0, -0.001, 0.0002],
            [-0.001, 0.0, 0.0002],
            [-0.001, 0.0002, 0.0],
        ]
    )

    eij = np.array(
        [
            [0.0, -0.00001, 0.00001],
            [-0.00001, 0.0, 0.00001],
            [-0.00001, 0.00001, 0.0],
        ]
    )

    rs = [0.92, 2.1055, 1.5]
    qs = [1.4, 1.972, 1.4]

    uniquac = UNIQUAC(qs, rs, aij, bij, cij, dij, eij)

    # =========================================================================
    # Models dictionary
    # =========================================================================
    models = {
        "NRTL": nrtl,
        "UNIFAC": unifac,
        "PSRK": psrk,
        "Dortmund": unifac_dortmund,
        "UNIQUAC": uniquac,
    }

    # =========================================================================
    # Test models
    # =========================================================================
    n = [15.9754, 3.125, 24.6721]
    temp = 320.0

    for line in data_lines:
        values = line.split(",")

        values = [v.strip() for v in values]

        model_name = values[0]
        model = models[model_name]

        thermoprops = [float(v) for v in values[1:]]

        ge, get, get2 = thermoprops[0:3]
        gen = thermoprops[3:6]
        getn = thermoprops[6:9]
        gen2 = np.reshape(thermoprops[9:18], (3, 3))

        (
            he,
            het,
        ) = thermoprops[18:20]
        hen = thermoprops[20:23]

        se, se_t = thermoprops[23:25]
        sen = thermoprops[25:28]

        lngamma = thermoprops[28:31]
        dlngamma_dt = thermoprops[31:34]
        dlngamma_dn = np.reshape(thermoprops[34:], (3, 3))

        # Test GE
        ge_v, derivatives = model.excess_gibbs(
            n, temp, dt=True, dt2=True, dn=True, dtn=True, dn2=True
        )

        ge_i = model.excess_gibbs(n, temp)

        assert np.isclose(ge_i, ge_v, rtol=1e-10)
        assert np.isclose(ge, ge_v, rtol=1e-10)
        assert np.isclose(get, derivatives["dt"], rtol=1e-10)
        assert np.isclose(get2, derivatives["dt2"], rtol=1e-10)
        assert np.allclose(gen, derivatives["dn"], rtol=1e-10)
        assert np.allclose(getn, derivatives["dtn"], rtol=1e-10)
        assert np.allclose(gen2, derivatives["dn2"], rtol=1e-10)

        # Test HE
        he_v, derivatives = model.excess_enthalpy(n, temp, dt=True, dn=True)

        he_i = model.excess_enthalpy(n, temp)

        assert np.isclose(he_i, he_v, rtol=1e-10)
        assert np.isclose(he, he_v, rtol=1e-10)
        assert np.isclose(het, derivatives["dt"], rtol=1e-10)
        assert np.allclose(hen, derivatives["dn"], rtol=1e-10)

        # Test SE
        se_v, derivatives = model.excess_entropy(n, temp, dt=True, dn=True)

        se_i = model.excess_entropy(n, temp)

        assert np.isclose(se_i, se_v, rtol=1e-10)
        assert np.isclose(se, se_v, rtol=1e-10)
        assert np.isclose(se_t, derivatives["dt"], rtol=1e-10)
        assert np.allclose(sen, derivatives["dn"], rtol=1e-10)

        # Test ln(gamma)
        lngamma_v, derivatives = model.ln_gamma(n, temp, dt=True, dn=True)

        lngamma_i = model.ln_gamma(n, temp)

        assert np.allclose(lngamma_i, lngamma_v, rtol=1e-10)
        assert np.allclose(lngamma, lngamma_v, rtol=1e-10)
        assert np.allclose(dlngamma_dt, derivatives["dt"], rtol=1e-10)
        assert np.allclose(dlngamma_dn, derivatives["dn"], rtol=1e-10)

    data_path.unlink()
