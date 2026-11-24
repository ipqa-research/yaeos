"""Tests for the constant volume translation Python interface."""

import gc

import numpy as np

import pytest

from yaeos import ConstantVolumeTranslation, PengRobinson76, QMR
from yaeos.constants import R
from yaeos.core import ArModel


def make_base():
    """Construct a binary methane/ethane model."""
    return PengRobinson76([190.564, 305.32], [45.99, 48.72], [0.01142, 0.0995])


def test_shift_and_thermodynamic_properties():
    """Compare corresponding translated and untranslated states."""
    base = make_base()
    c = np.array([0.012, 0.021])
    model = ConstantVolumeTranslation(base, c)
    n = np.array([0.8, 1.5])
    temperature = 330.0
    volume = 0.8
    base_volume = volume + n @ c

    assert isinstance(model, ArModel)
    assert model.size() == 2
    pressure = base.pressure(n, base_volume, temperature)
    assert model.pressure(n, volume, temperature) == pytest.approx(pressure)

    # The ideal-reference correction is essential in the residual energy.
    base_ar = base.helmholtz_residual_vt(n, base_volume, temperature)
    ar = model.helmholtz_residual_vt(n, volume, temperature)
    expected = base_ar + n.sum() * R * temperature * np.log(
        volume / base_volume
    )
    assert ar == pytest.approx(expected, rel=1e-12, abs=1e-12)

    for root in ["liquid", "vapor", "stable"]:
        base_v = base.volume(n, 1.0, 200.0, root=root)
        shifted_v = model.volume(n, 1.0, 200.0, root=root)
        assert shifted_v > 0
        assert shifted_v == pytest.approx(base_v - n @ c, rel=1e-12)
        assert model.pressure(n, shifted_v, 200.0) == pytest.approx(
            1.0, rel=1e-6
        )


def test_helmholtz_optional_derivatives():
    """Individual derivative buffers match the combined call."""
    model = ConstantVolumeTranslation(make_base(), [0.012, 0.021])
    n = [0.8, 1.5]
    options = ["dt", "dv", "dn", "dtv", "dv2", "dt2", "dvn", "dtn", "dn2"]
    energy, derivatives = model.helmholtz_residual_vt(
        n, 0.8, 330.0, **{name: True for name in options}
    )
    assert np.isfinite(energy)

    for name in options:
        individual_energy, individual = model.helmholtz_residual_vt(
            n, 0.8, 330.0, **{name: True}
        )
        assert individual_energy == pytest.approx(energy)
        np.testing.assert_allclose(individual[name], derivatives[name])
        for other in options:
            if other != name:
                assert individual[other] is None

    _, combined = model.helmholtz_residual_vt(
        n, 0.8, 330.0, dtv=True, dvn=True
    )
    np.testing.assert_allclose(combined["dtv"], derivatives["dtv"])
    np.testing.assert_allclose(combined["dvn"], derivatives["dvn"])
    assert derivatives["dn"].shape == (2,)
    assert derivatives["dn2"].shape == (2, 2)


def test_copy_and_registry_lifetime():
    """A translated model owns its base, independently of Python handles."""
    base = make_base()
    shifts = np.array([0.012, 0.021])
    model = ConstantVolumeTranslation(base, shifts)
    expected = model.pressure([0.8, 1.5], 0.8, 330.0)

    shifts[:] = 0.0
    base.set_mixrule(
        QMR(kij=[[0.0, 0.4], [0.4, 0.0]], lij=[[0.0, 0.0], [0.0, 0.0]])
    )
    assert model.pressure([0.8, 1.5], 0.8, 330.0) == pytest.approx(expected)

    del base
    gc.collect()
    replacement = make_base()
    assert model.pressure([0.8, 1.5], 0.8, 330.0) == pytest.approx(expected)
    del replacement

    for _ in range(10):
        temporary = ConstantVolumeTranslation(make_base(), [0.01, 0.02])
        assert temporary.pressure([0.8, 1.5], 0.8, 330.0) > 0
        del temporary

    assert model.pressure([0.8, 1.5], 0.8, 330.0) == pytest.approx(expected)


@pytest.mark.parametrize(
    "c", [[], [0.01], [0.01, 0.02, 0.03], [[0.01, 0.02]], [0.01, np.nan]]
)
def test_invalid_shifts(c):
    """Reject invalid coefficients before calling the native wrapper."""
    with pytest.raises(ValueError):
        ConstantVolumeTranslation(make_base(), c)


def test_invalid_base():
    """Only residual Helmholtz models can be translated."""
    with pytest.raises(TypeError):
        ConstantVolumeTranslation(object(), [0.01])


def test_fortran_representation_declarations_and_setup():
    """Declare models before initialization, including nested shifts."""
    base = make_base()
    base.set_mixrule(
        QMR(kij=[[0.0, 0.1], [0.1, 0.0]], lij=[[0.0, 0.0], [0.0, 0.0]])
    )
    inner = ConstantVolumeTranslation(base, [0.01, 0.02])
    model = ConstantVolumeTranslation(inner, [0.003, 0.004])
    declarations = model._model_params_declaration_as_str()
    setup = model._model_params_as_str()

    assert declarations.count("parameter :: nc=") == 1
    assert "type(CubicEoS) :: base_model" in declarations
    assert "type(ConstantVolumeTranslation) :: base_model_base" in declarations
    assert "type(ConstantVolumeTranslation) :: ar_model" in declarations
    assert "block" not in declarations + setup
    assert "call base_model%set_mixrule(mixrule)" in setup
    assert "base_model_base = ConstantVolumeTranslation(base_model," in setup
    assert "ar_model = ConstantVolumeTranslation(base_model_base," in setup


@pytest.mark.parametrize(
    "method, expected",
    [
        ("psrk_correlation", [0.0384890688, -0.0022789312]),
        ("vtpr_correlation", [-0.02354688, 0.01538208]),
    ],
)
def test_correlations(method, expected, monkeypatch):
    """Evaluate both correlations without constructing a model."""

    def fail_init(*args, **kwargs):
        raise AssertionError("correlation must not construct a model")

    monkeypatch.setattr(ConstantVolumeTranslation, "__init__", fail_init)
    # RTc/Pc = 1 L/mol; Zc is 0.2 and 0.3, respectively.
    result = getattr(ConstantVolumeTranslation, method)(
        [300.0, 300.0], [300.0 * R, 300.0 * R], [0.2, 0.3]
    )
    assert isinstance(result, np.ndarray)
    np.testing.assert_allclose(result, expected, atol=1e-12, rtol=0)


@pytest.mark.parametrize("method", ["psrk_correlation", "vtpr_correlation"])
@pytest.mark.parametrize(
    "tc, pc, vc",
    [
        ([], [], []),
        ([300.0], [40.0, 50.0], [0.1]),
        ([[300.0]], [40.0], [0.1]),
        ([0.0], [40.0], [0.1]),
        ([300.0], [-1.0], [0.1]),
        ([300.0], [40.0], [np.nan]),
    ],
)
def test_correlations_reject_invalid_properties(method, tc, pc, vc):
    """Reject invalid or incompatible pure-component properties."""
    with pytest.raises(ValueError):
        getattr(ConstantVolumeTranslation, method)(tc, pc, vc)
