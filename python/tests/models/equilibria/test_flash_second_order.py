"""Tests for the second-order (Michelsen) two-phase flash.

The second-order flash minimizes directly the Gibbs energy (fixed ``T`` and
``P``) or the Helmholtz energy (fixed ``T``) instead of solving the
successive-substitution fixed point.  These tests check that the Python
bindings reach the same split as the classic successive-substitution flash
and the reference values already used by the Fortran test suite.
"""

import numpy as np

import pytest

from yaeos import (
    NRTL,
    UNIQUAC,
    UNIFACVLE,
    PengRobinson76,
    PengRobinson78,
    QMR,
    SoaveRedlichKwong,
)

# Methane/n-butane, used by the Fortran ``test_flash_pt``/``test_flash_pt_2``.
TC = [190.564, 425.12]
PC = [45.99, 37.96]
W = [0.0115478, 0.200164]


def _assert_same_split(result, reference, atol=1e-5):
    """Compare a second-order flash against the successive substitutions."""
    assert result["beta"] > 0.0
    np.testing.assert_allclose(result["x"], reference["x"], atol=atol, rtol=0)
    np.testing.assert_allclose(result["y"], reference["y"], atol=atol, rtol=0)
    assert np.isclose(result["beta"], reference["beta"], atol=atol, rtol=0)


@pytest.mark.parametrize(
    "model_class",
    [PengRobinson76, PengRobinson78, SoaveRedlichKwong],
)
@pytest.mark.parametrize("pressure", [20.0, 40.0, 60.0])
def test_ar_models_match_successive_substitution(model_class, pressure):
    """``flash_pt_second_order`` must reproduce ``flash_pt``."""
    model = model_class(TC, PC, W)
    z = [0.4, 0.6]
    temperature = 294.0

    reference = model.flash_pt(z, pressure, temperature)
    result = model.flash_pt_second_order(z, pressure, temperature)

    _assert_same_split(result, reference)

    # Volumes are reported as molar volumes of each phase.
    assert np.isclose(result["Vx"], reference["Vx"], atol=1e-5, rtol=0)
    assert np.isclose(result["Vy"], reference["Vy"], atol=1e-5, rtol=0)
    assert result["P"] == pressure
    assert result["T"] == temperature


def test_ar_wilson_and_k0_initializations_agree():
    """Both default (Wilson) and explicit ``k0`` must reach the same split."""
    model = PengRobinson76(TC, PC, W)
    z = np.array([0.4, 0.6])
    pressure = 60.0
    temperature = 294.0

    reference = model.flash_pt(z, pressure, temperature)

    k0 = (np.array(PC) / pressure) * np.exp(
        5.373 * (1.0 + np.array(W)) * (1.0 - np.array(TC) / temperature)
    )
    from_k0 = model.flash_pt_second_order(z, pressure, temperature, k0=k0)
    from_wilson = model.flash_pt_second_order(z, pressure, temperature)

    _assert_same_split(from_k0, reference)
    _assert_same_split(from_wilson, reference)


def test_ar_ternary_pr78():
    """Ternary Peng-Robinson 78 case from the Fortran test suite."""
    tc = [190.564, 304.1282, 768.0]
    pc = [45.992, 73.773, 10.7]
    w = [0.01142, 0.22394, 0.8805]

    kij = np.zeros((3, 3))
    kij[0, 1] = kij[1, 0] = 0.1
    kij[1, 2] = kij[2, 1] = 0.2

    model = PengRobinson78(tc, pc, w, QMR(kij, np.zeros((3, 3))))

    z = np.array([0.2, 0.4, 0.4])
    z = z / np.sum(z)
    temperature = 260.0
    pressure = 150.0

    reference = model.flash_pt(z, pressure, temperature)
    result = model.flash_pt_second_order(z, pressure, temperature)

    _assert_same_split(result, reference)
    assert np.isclose(result["Vx"], reference["Vx"], atol=1e-5, rtol=0)
    assert np.isclose(result["Vy"], reference["Vy"], atol=1e-5, rtol=0)


def test_ar_liquid_light_root():
    """A liquid-liquid split with the light phase forced to the liquid root."""
    tc = [500.0, 500.0]
    pc = [30.0, 30.0]
    w = [0.3, 0.3]

    kij = np.zeros((2, 2))
    kij[0, 1] = kij[1, 0] = 0.2
    model = PengRobinson76(tc, pc, w, QMR(kij, np.zeros((2, 2))))

    z = np.array([0.5, 0.5])
    temperature = 300.0
    pressure = 1.0

    stability, _ = model.stability_analysis(z, pressure, temperature)
    k0 = stability["w"] / z

    reference = model.flash_pt(z, pressure, temperature, k0=k0)
    result = model.flash_pt_second_order(
        z, pressure, temperature, k0=k0, light_root="liquid"
    )

    _assert_same_split(result, reference)

    # Both phases are liquid-like.
    assert result["Vx"] < 0.2
    assert result["Vy"] < 0.2


def test_ar_non_converged_returns_invalid_beta():
    """A non-converged second-order flash reports ``beta == -1``."""
    model = PengRobinson76(TC, PC, W)
    z = [0.4, 0.6]

    result = model.flash_pt_second_order(z, 60.0, 294.0, max_iters=0)

    assert result["beta"] == -1.0
    np.testing.assert_allclose(result["x"], z)
    np.testing.assert_allclose(result["y"], z)


def test_ge_unifac_matches_flash_t():
    """UNIFAC liquid-liquid split from the Fortran test suite."""
    model = UNIFACVLE([{1: 1, 42: 1}, {1: 2, 2: 4}, {16: 1}])

    z = np.array([0.2, 0.5, 0.3])
    temperature = 250.0

    k0 = np.array([10.699186891095467, 5.685729557771e-3, 476.8624745978])

    reference = model.flash_t(z, temperature, k0=k0)
    from_k0 = model.flash_t_second_order(z, temperature, k0=k0)

    # The default initialization uses a stability analysis.
    from_stability = model.flash_t_second_order(z, temperature)

    _assert_same_split(from_k0, reference)
    _assert_same_split(from_stability, reference)


def test_ge_nrtl_matches_flash_t():
    """Strongly non-ideal binary NRTL mixture with two liquid phases."""
    a = np.zeros((2, 2))
    b = np.zeros((2, 2))
    c = np.zeros((2, 2))
    a[0, 1] = a[1, 0] = 5.0
    c[0, 1] = c[1, 0] = 0.3

    model = NRTL(a, b, c)

    z = np.array([0.3, 0.7])
    temperature = 300.0
    k0 = np.array([405.70535832398622, 2.4648429478662474e-3])

    reference = model.flash_t(z, temperature, k0=k0)
    result = model.flash_t_second_order(z, temperature, k0=k0)

    _assert_same_split(result, reference)


def test_ge_uniquac_successive_substitution_failure_point():
    """A point where the successive-substitution flash does not converge.

    The reference values come from the Gibbs-energy minimization and are also
    checked in the Fortran test suite.
    """
    nc = 7

    # Component order: F, H2O2, PF, H2O, oil, ESBO, polyol.
    n = np.array(
        [
            8.69085708798058754e00,
            1.04599576389094473e00,
            8.06134619967180610e-01,
            2.12080199832480858e01,
            4.53403810029077592e-01,
            4.54221490480090750e-01,
            8.39249341818703722e-14,
        ]
    )
    z = n / np.sum(n)
    temperature = 333.15

    qs = np.array([1.532, 2.4, 2.388, 1.4, 31.347, 30.492, 41.232])
    rs = np.array([1.528, 2.0, 2.242, 0.92, 38.7157, 38.8207, 48.8112])

    aij = np.zeros((nc, nc))
    aij[0, 3] = 0.0650752585
    aij[3, 0] = -1.82683869

    bij = np.zeros((nc, nc))
    bij[0, :] = [
        0.0,
        151.0,
        135.49896199,
        59.3246124,
        -11095.62551048,
        -3246.2424174,
        -995.317802,
    ]
    bij[1, :] = [
        -199.0,
        0.0,
        -101.29550998,
        -11.76497438,
        -62.99691978,
        -1940.97729389,
        -675.950034,
    ]
    bij[2, :] = [
        52.89079805,
        -3.03767903,
        0.0,
        107.87953158,
        -1444.81779555,
        -1113.00536155,
        -787.3298,
    ]
    bij[3, :] = [
        947.601438,
        173.79237634,
        -210.13173272,
        0.0,
        158.23356021,
        190.74747368,
        180.790188,
    ]
    bij[4, :] = [
        21.69729307,
        -169.18769991,
        207.5615499,
        -1252.35267809,
        0.0,
        0.0,
        0.0,
    ]
    bij[5, :] = [
        346.55860276,
        142.25776765,
        414.04790322,
        -1308.66514598,
        0.0,
        0.0,
        0.0,
    ]
    bij[6, :] = [
        264.640828,
        37.832309,
        340.930152,
        -1362.762912,
        0.0,
        0.0,
        0.0,
    ]

    zeros = np.zeros((nc, nc))
    model = UNIQUAC(qs, rs, aij, bij, zeros, zeros, zeros)

    k0 = np.array(
        [
            1.04947543188346137e00,
            5.59102656162917544e-01,
            2.05883043639840357e00,
            2.63674340469321811e-01,
            4.31976399332579896e07,
            1.79109736675408584e01,
            4.05540135724048734e-02,
        ]
    )

    x_ref = np.array([0.2656, 0.0326, 0.0236, 0.6700])
    y_ref = np.array([0.2787, 0.0182, 0.0487, 0.1767])
    beta_ref = 0.0418376

    result = model.flash_t_second_order(z, temperature, k0=k0)

    assert result["beta"] > 0.0
    assert np.isclose(result["beta"], beta_ref, atol=1e-3, rtol=0)
    np.testing.assert_allclose(result["x"][:4], x_ref, atol=1e-3, rtol=0)
    np.testing.assert_allclose(result["y"][:4], y_ref, atol=1e-3, rtol=0)
    assert np.max(np.abs(result["x"] - result["y"])) > 1e-2

    # Equal activities are the equilibrium condition for an activity model.
    activity_x = np.exp(model.ln_gamma(result["x"], temperature)) * result["x"]
    activity_y = np.exp(model.ln_gamma(result["y"], temperature)) * result["y"]
    np.testing.assert_allclose(activity_x, activity_y, atol=1e-6, rtol=0)

    # The stability-analysis initialization reaches the same split.
    from_stability = model.flash_t_second_order(z, temperature)

    np.testing.assert_allclose(
        from_stability["x"][:4], x_ref, atol=1e-3, rtol=0
    )
    np.testing.assert_allclose(
        from_stability["y"][:4], y_ref, atol=1e-3, rtol=0
    )
    assert np.isclose(from_stability["beta"], beta_ref, atol=1e-3, rtol=0)


def test_ge_trivial_solution_rejected():
    """``k0 = 1`` is the trivial solution and must be rejected."""
    model = UNIFACVLE([{1: 1, 42: 1}, {1: 2, 2: 4}, {16: 1}])
    z = np.array([0.2, 0.5, 0.3])

    result = model.flash_t_second_order(z, 250.0, k0=np.ones_like(z))

    assert result["beta"] == -1.0
    np.testing.assert_allclose(result["x"], z)
    np.testing.assert_allclose(result["y"], z)
