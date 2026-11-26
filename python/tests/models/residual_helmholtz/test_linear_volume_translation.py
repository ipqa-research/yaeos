"""Tests for the linear volume translation Python interface."""

import gc

import numpy as np

import pytest

from yaeos import (
    ConstantVolumeTranslation,
    LinearVolumeTranslation,
    PengRobinson76,
    QMR,
)
from yaeos.constants import R
from yaeos.core import ArModel


def make_base():
    """Construct a binary methane/ethane model."""
    return PengRobinson76([190.564, 305.32], [45.99, 48.72], [0.01142, 0.0995])


@pytest.mark.parametrize("temperature", [200.0, 330.0, 450.0])
def test_corresponding_states(temperature):
    """Check energy, pressure and vapor volume at corresponding states."""
    base = make_base()
    a = np.array([0.012, 0.021])
    b = np.array([4e-5, -2e-5])
    model = LinearVolumeTranslation(base, a, b)
    n = np.array([0.8, 1.5])
    volume = 0.8
    shift = n @ (a + b * temperature)
    base_volume = volume + shift

    assert isinstance(model, ArModel)
    assert model.size() == 2
    assert model.pressure(n, volume, temperature) == pytest.approx(
        base.pressure(n, base_volume, temperature), rel=1e-12
    )
    base_ar = base.helmholtz_residual_vt(n, base_volume, temperature)
    ar = model.helmholtz_residual_vt(n, volume, temperature)
    expected = base_ar + n.sum() * R * temperature * np.log(
        volume / base_volume
    )
    assert ar == pytest.approx(expected, rel=1e-12, abs=1e-12)

    base_v = base.volume(n, 1.0, temperature, root="vapor")
    shifted_v = model.volume(n, 1.0, temperature, root="vapor")
    assert base_v > 0 and shifted_v > 0
    assert shifted_v == pytest.approx(base_v - shift, rel=1e-12)
    assert model.pressure(n, shifted_v, temperature) == pytest.approx(
        1.0, rel=1e-6
    )


def test_optional_helmholtz_derivatives():
    """Individual optional outputs agree with the combined call."""
    model = LinearVolumeTranslation(make_base(), [0.012, 0.021], [4e-5, -2e-5])
    names = ["dt", "dv", "dn", "dtv", "dv2", "dt2", "dvn", "dtn", "dn2"]
    energy, all_derivatives = model.helmholtz_residual_vt(
        [0.8, 1.5], 0.8, 330.0, **{name: True for name in names}
    )
    for name in names:
        value, individual = model.helmholtz_residual_vt(
            [0.8, 1.5], 0.8, 330.0, **{name: True}
        )
        assert value == pytest.approx(energy)
        np.testing.assert_allclose(individual[name], all_derivatives[name])
        assert all(
            individual[other] is None for other in names if other != name
        )
    assert all_derivatives["dn"].shape == (2,)
    assert all_derivatives["dn2"].shape == (2, 2)


def test_temperature_derivatives_against_differences():
    """Exercise the nonzero temperature dependence through the binding."""
    model = LinearVolumeTranslation(make_base(), [0.012, 0.021], [4e-5, -2e-5])
    n = [0.8, 1.5]
    v = 0.8
    t = 330.0
    h = 0.01
    _, derivatives = model.helmholtz_residual_vt(n, v, t, dt=True, dt2=True)
    plus = model.helmholtz_residual_vt(n, v, t + h)
    minus = model.helmholtz_residual_vt(n, v, t - h)
    _, dplus = model.helmholtz_residual_vt(n, v, t + h, dt=True)
    _, dminus = model.helmholtz_residual_vt(n, v, t - h, dt=True)
    assert derivatives["dt"] == pytest.approx(
        (plus - minus) / (2 * h), rel=1e-7
    )
    assert derivatives["dt2"] == pytest.approx(
        (dplus["dt"] - dminus["dt"]) / (2 * h), rel=1e-7
    )


def test_zero_slopes_match_constant():
    """Zero slopes reproduce all derivatives of the constant translation."""
    base = make_base()
    linear = LinearVolumeTranslation(base, [0.012, 0.021], [0.0, 0.0])
    constant = ConstantVolumeTranslation(base, [0.012, 0.021])
    names = ["dt", "dv", "dn", "dtv", "dv2", "dt2", "dvn", "dtn", "dn2"]
    opts = {name: True for name in names}
    lvalue, lder = linear.helmholtz_residual_vt([0.8, 1.5], 0.8, 330, **opts)
    cvalue, cder = constant.helmholtz_residual_vt([0.8, 1.5], 0.8, 330, **opts)
    assert lvalue == pytest.approx(cvalue, rel=1e-12)
    for name in names:
        np.testing.assert_allclose(
            lder[name], cder[name], rtol=1e-12, atol=1e-12
        )


def test_copy_and_registry_lifetime():
    """Neither mutating inputs nor releasing the base changes the copy."""
    base = make_base()
    a = np.array([0.012, 0.021])
    b = np.array([4e-5, -2e-5])
    model = LinearVolumeTranslation(base, a, b)
    expected = model.pressure([0.8, 1.5], 0.8, 330.0)
    a[:] = 0
    b[:] = 0
    base.set_mixrule(
        QMR(kij=[[0.0, 0.4], [0.4, 0.0]], lij=[[0.0, 0.0], [0.0, 0.0]])
    )
    del base
    gc.collect()
    for _ in range(10):
        temporary = LinearVolumeTranslation(make_base(), [0.01, 0.02], [0, 0])
        assert temporary.pressure([0.8, 1.5], 0.8, 330.0) > 0
        del temporary
    assert model.pressure([0.8, 1.5], 0.8, 330.0) == pytest.approx(expected)
    assert not model.a.flags.writeable and not model.b.flags.writeable


@pytest.mark.parametrize("name", ["a", "b"])
@pytest.mark.parametrize(
    "invalid", [[], [0.1], [0, 0, 0], [[0, 0]], [0, np.nan], [0, np.inf]]
)
def test_invalid_coefficients(name, invalid):
    """Reject invalid vectors before entering the native code."""
    args = {"a": [0.01, 0.02], "b": [1e-5, -2e-5]}
    args[name] = invalid
    with pytest.raises(ValueError):
        LinearVolumeTranslation(make_base(), **args)


def test_invalid_base():
    """Only residual Helmholtz models can be translated."""
    with pytest.raises(TypeError):
        LinearVolumeTranslation(object(), [0.01], [0.001])


def test_nested_fortran_representation():
    """Declare both models before initializing, including nested shifts."""
    base = ConstantVolumeTranslation(make_base(), [0.01, 0.02])
    model = LinearVolumeTranslation(base, [0.001, 0.002], [1e-5, -2e-5])
    declarations = model._model_params_declaration_as_str()
    setup = model._model_params_as_str()
    assert declarations.count("parameter :: nc=") == 1
    assert "type(CubicEoS) :: base_model" in declarations
    assert "type(ConstantVolumeTranslation) :: base_model_base" in declarations
    assert "type(LinearVolumeTranslation) :: ar_model" in declarations
    assert "block" not in declarations + setup
    assert "base_model_base = ConstantVolumeTranslation(base_model," in setup
    assert "ar_model = LinearVolumeTranslation(base_model_base," in setup
