"""Volume translations for residual Helmholtz models."""

import re

import numpy as np

from yaeos.constants import R
from yaeos.core import ArModel
from yaeos.lib import yaeos_c


def _critical_properties(
    critical_temperatures, critical_pressures, critical_volumes
):
    """Validate and convert pure-component critical properties."""
    arrays = tuple(
        np.asarray(values, dtype=float)
        for values in (
            critical_temperatures,
            critical_pressures,
            critical_volumes,
        )
    )
    if any(a.ndim != 1 or a.size == 0 for a in arrays):
        raise ValueError("critical properties must be nonempty vectors")
    if any(a.shape != arrays[0].shape for a in arrays):
        raise ValueError("critical property vectors must have the same length")
    if any(not np.all(np.isfinite(a) & (a > 0)) for a in arrays):
        raise ValueError("critical properties must be positive and finite")
    return arrays


class ConstantVolumeTranslation(ArModel):
    r"""Apply constant component volume shifts to a residual Helmholtz model.

    The extensive shift is :math:`C=\sum_i n_i c_i`. The base EOS is
    evaluated at :math:`W=V+C`, where ``V`` is the physical total volume.
    Positive shifts therefore reduce the volume at fixed pressure and
    temperature. The base EOS is copied when this model is constructed;
    subsequent changes to the original model do not change the copy.

    Parameters
    ----------
    base : ArModel
        Initialized residual Helmholtz model to translate.
    c : array_like
        One constant molar shift per component, in base-model order [L/mol].

    Examples
    --------
    >>> from yaeos import PengRobinson76, ConstantVolumeTranslation
    >>> base = PengRobinson76([190.564], [45.99], [0.01142])
    >>> model = ConstantVolumeTranslation(base, [0.01])
    >>> volume = model.volume([1.0], 1.0, 300.0, root="vapor")
    """

    def __init__(self, base: ArModel, c):
        if not isinstance(base, ArModel):
            raise TypeError("base must be an ArModel")

        coefficients = np.asarray(c, dtype=np.float64)
        if coefficients.ndim != 1 or coefficients.size != base.size():
            raise ValueError("c must contain one shift per base component")
        if not np.all(np.isfinite(coefficients)):
            raise ValueError("c must contain finite shifts")

        self.c = coefficients.copy()
        self.c.flags.writeable = False
        self._nc = base.size()
        # Keep the comparison-test representation consistent with the copy
        # taken by Fortran, even if the original model is later reconfigured.
        self._base_declaration = base._model_params_declaration_as_str()
        self._base_setup = base._model_params_as_str()
        self._base_name = "base_model"
        while re.search(
            rf"\b{self._base_name}\b", self._base_declaration
        ) or re.search(rf"\b{self._base_name}\b", self._base_setup):
            self._base_name += "_base"
        self.id = yaeos_c.constant_volume_translation(base.id, self.c)

    def size(self) -> int:
        """Return the number of components."""
        return self._nc

    @classmethod
    def psrk_correlation(
        cls, critical_temperatures, critical_pressures, critical_volumes
    ):
        r"""Calculate SRK shifts using the PSRK correlation.

        Horstmann et al. (2005), Fluid Phase Equilibria 227, 157–164,
        equations 8–9. The experimental critical compressibility is
        :math:`Z_c=P_c v_c/(RT_c)`, and
        :math:`c=0.40768(RT_c/P_c)(0.29441-Z_c)`.

        Parameters
        ----------
        critical_temperatures : array_like
            Pure-component critical temperatures [K].
        critical_pressures : array_like
            Pure-component critical pressures [bar].
        critical_volumes : array_like
            Experimental critical molar volumes [L/mol].

        Returns
        -------
        numpy.ndarray
            Component shifts [L/mol], in input order.
        """
        tc, pc, vc = _critical_properties(
            critical_temperatures, critical_pressures, critical_volumes
        )
        zc = pc * vc / (R * tc)
        return 0.40768 * R * tc / pc * (0.29441 - zc)

    @classmethod
    def vtpr_correlation(
        cls, critical_temperatures, critical_pressures, critical_volumes
    ):
        r"""Calculate PR76 shifts using the Ahlers–Gmehling VTPR correlation.

        Ahlers and Gmehling (2001), Fluid Phase Equilibria 191, 177–188,
        equations 4 and 12–13. The experimental critical compressibility is
        :math:`Z_c=P_c v_c/(RT_c)`, and
        :math:`c=0.252(RT_c/P_c)(1.5448Z_c-0.4024)`.
        These shifts alone do not construct a VTPR model.

        Parameters
        ----------
        critical_temperatures : array_like
            Pure-component critical temperatures [K].
        critical_pressures : array_like
            Pure-component critical pressures [bar].
        critical_volumes : array_like
            Experimental critical molar volumes [L/mol].

        Returns
        -------
        numpy.ndarray
            Component shifts [L/mol], in input order.
        """
        tc, pc, vc = _critical_properties(
            critical_temperatures, critical_pressures, critical_volumes
        )
        zc = pc * vc / (R * tc)
        return 0.252 * R * tc / pc * (1.5448 * zc - 0.4024)

    def _model_params_declaration_as_str(self) -> str:
        """Return Fortran declarations for the translated model."""
        declaration = re.sub(
            r"\bar_model\b", self._base_name, self._base_declaration
        )
        return declaration + "\ntype(ConstantVolumeTranslation) :: ar_model\n"

    def _model_params_as_str(self) -> str:
        """Return Fortran setup code for the copied EOS and its translation."""
        setup = re.sub(r"\bar_model\b", self._base_name, self._base_setup)
        shifts = ", ".join(f"{value}_pr" for value in self.c)
        return (
            f"{setup}\n"
            f"ar_model = ConstantVolumeTranslation({self._base_name}, "
            f"[{shifts}])\n"
        )

    def __del__(self):
        """Release the registry slot if construction succeeded."""
        if hasattr(self, "id"):
            yaeos_c.make_available_ar_models_list(self.id)


class LinearVolumeTranslation(ArModel):
    r"""Apply linear temperature-dependent component volume shifts.

    Each shift is :math:`c_i(T)=a_i+b_iT`, and the base EOS is evaluated
    at :math:`W=V+\sum_i n_i c_i(T)`. Positive shifts reduce the physical
    volume at fixed pressure and temperature. The configured base model is
    copied on construction; later changes to the original do not affect it.

    Parameters
    ----------
    base : ArModel
        Initialized residual Helmholtz model to translate.
    a : array_like
        Component shift intercepts in base-model order [L/mol].
    b : array_like
        Component shift slopes in base-model order [L/(mol K)].

    Examples
    --------
    >>> from yaeos import PengRobinson76, LinearVolumeTranslation
    >>> base = PengRobinson76([190.564], [45.99], [0.01142])
    >>> model = LinearVolumeTranslation(base, [0.01], [-1e-5])
    >>> volume = model.volume([1.0], 1.0, 300.0, root="vapor")
    """

    def __init__(self, base: ArModel, a, b):
        if not isinstance(base, ArModel):
            raise TypeError("base must be an ArModel")

        self._nc = base.size()
        for name, values in (("a", a), ("b", b)):
            coefficients = np.asarray(values, dtype=np.float64)
            if coefficients.ndim != 1 or coefficients.size != self._nc:
                raise ValueError(
                    f"{name} must contain one coefficient per base component"
                )
            if not np.all(np.isfinite(coefficients)):
                raise ValueError(f"{name} must contain finite coefficients")
            coefficients = coefficients.copy()
            coefficients.flags.writeable = False
            setattr(self, name, coefficients)

        # Snapshot the base representation as well as the native model.
        self._base_declaration = base._model_params_declaration_as_str()
        self._base_setup = base._model_params_as_str()
        self._base_name = "base_model"
        while re.search(
            rf"\b{self._base_name}\b", self._base_declaration
        ) or re.search(rf"\b{self._base_name}\b", self._base_setup):
            self._base_name += "_base"
        self.id = yaeos_c.linear_volume_translation(base.id, self.a, self.b)

    def size(self) -> int:
        """Return the number of components."""
        return self._nc

    def _model_params_declaration_as_str(self) -> str:
        """Return Fortran declarations for both models."""
        declaration = re.sub(
            r"\bar_model\b", self._base_name, self._base_declaration
        )
        return declaration + "\ntype(LinearVolumeTranslation) :: ar_model\n"

    def _model_params_as_str(self) -> str:
        """Return Fortran setup code for both models."""
        setup = re.sub(r"\bar_model\b", self._base_name, self._base_setup)
        a = ", ".join(f"{value}_pr" for value in self.a)
        b = ", ".join(f"{value}_pr" for value in self.b)
        return (
            f"{setup}\n"
            f"ar_model = LinearVolumeTranslation({self._base_name}, "
            f"[{a}], [{b}])\n"
        )

    def __del__(self):
        """Release the registry slot if construction succeeded."""
        if hasattr(self, "id"):
            yaeos_c.make_available_ar_models_list(self.id)
