<div align="center">

[![Fortran](https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat)](https://fortran-lang.org)
[![fpm](https://img.shields.io/badge/fpm-Fortran_package_manager-734f96)](https://fpm.fortran-lang.org)
[![FORD Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://ipqa-research.github.io/yaeos/)
[![Python Docs](https://img.shields.io/badge/Python-Docs-3776AB?logo=python&logoColor=white)](https://ipqa-research.github.io/yaeos/page/python-api/index.html)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://github.com/ipqa-research/yaeos/blob/main/LICENSE)
[![CI](https://github.com/ipqa-research/yaeos/actions/workflows/CI.yml/badge.svg)](https://github.com/ipqa-research/yaeos/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ipqa-research/yaeos/graph/badge.svg?token=IDJYKV8XK6)](https://codecov.io/gh/ipqa-research/yaeos)

<p align="center">
<a target="_blank" href="https://colab.research.google.com/github/ipqa-research/yaeos/blob/main/notebooks/run_online.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a> <br>
<img src="media/logo.png" width="50%" alt="yaeos logo"/>
</p>

</div>

**yaeos** (Yet Another Equation Of State library) computes thermodynamic
properties and phase equilibria with equations of state, in modern Fortran and
Python. It focuses on robust algorithms for phase diagrams and on an easy
framework for adding new models, with analytical derivatives or automatic
differentiation.

- **Phase diagrams that close.** Continuation methods trace PT, Px and Tx
  envelopes through critical points, including three-phase and multiphase
  boundaries, critical lines and binary global phase diagrams.
- **Models you can extend.** Write the residual Helmholtz or excess Gibbs
  energy and get fugacities, residual properties and every equilibrium routine
  for free.
- **Fortran speed, Python reach.** A readable Fortran core built with
  [fpm](https://fpm.fortran-lang.org), wrapped as a `pip install yaeos`
  package for notebooks, parameter fitting and plotting.

## Documentation

- [User guide and Fortran API](https://ipqa-research.github.io/yaeos/)
- [Python API](https://ipqa-research.github.io/yaeos/page/python-api/index.html)
- [Adding new models](https://ipqa-research.github.io/yaeos/page/usage/newmodels/index.html)
- [Developer handbook](https://ipqa-research.github.io/yaeos/page/developer/index.html)

## Installation

### Fortran

yaeos is an [fpm](https://fpm.fortran-lang.org) package. It needs BLAS and
LAPACK; on Debian/Ubuntu:

```bash
sudo apt install libblas-dev liblapack-dev
```

Then add it as a dependency in your project's `fpm.toml`:

```toml
[dependencies]
yaeos = {git="https://github.com/ipqa-research/yaeos"}
```

### Python

Binary wheels are published on PyPI:

```bash
pip install yaeos
```

## Quick start

Build a Peng-Robinson model for methane + n-butane, then evaluate a pressure,
a bubble point and the PT phase envelope. yaeos uses K, bar, L and mol
throughout.

### Fortran

```fortran
program quickstart
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, PTEnvel2, &
      saturation_pressure, pt_envelope_2ph
   implicit none
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: bubble
   type(PTEnvel2) :: envelope
   real(pr) :: z(2)
   real(pr) :: P

   ! Methane + n-butane with Peng-Robinson (1976): Tc [K], Pc [bar], omega
   model = PengRobinson76(tc=[190.564_pr, 425.12_pr], &
      pc=[45.99_pr, 37.96_pr], w=[0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]  ! mol

   call model%pressure(z, V=1.0_pr, T=300.0_pr, P=P)
   bubble = saturation_pressure(model, z, T=300.0_pr, kind="bubble")
   envelope = pt_envelope_2ph(model, z, &
      saturation_pressure(model, z, T=200.0_pr, kind="dew"))

   print *, "P(1 L, 300 K) =", P, "bar"
   print *, "Bubble P(300 K) =", bubble%P, "bar"
   print *, "Critical point:", envelope%cps(1)%T, "K", envelope%cps(1)%P, "bar"
end program quickstart
```

### Python

```python
import numpy as np
from yaeos import PengRobinson76

# Methane + n-butane with Peng-Robinson (1976)
model = PengRobinson76(
    critical_temperatures=[190.564, 425.12],  # K
    critical_pressures=[45.99, 37.96],  # bar
    acentric_factors=[0.0115478, 0.200164],
)
z = np.array([0.3, 0.7])  # mol

P = model.pressure(z, volume=1.0, temperature=300.0)
bubble = model.saturation_pressure(z, temperature=300.0, kind="bubble")
envelope = model.phase_envelope_pt(z, kind="dew", t0=200.0, p0=1.0)

print(f"P(1 L, 300 K) = {P:.4f} bar")
print(f"Bubble P(300 K) = {bubble['P']:.4f} bar")
print("Critical point:", envelope.critical_temperatures, envelope.critical_pressures)
envelope.plot()  # requires matplotlib
```

Both print a pressure of 15.91 bar, a bubble pressure of 57.37 bar and a
critical point near 401.5 K and 68.1 bar.

### Examples

More Fortran examples live in [`example/`](example). From a clone of this
repository, run a registered example by name:

```bash
fpm run --example basics
```

Available names (see `fpm.toml`): `basics`, `saturation_points`,
`phase_split`, `phase_envelope`, `new_alpha_function`, `pure_psat`,
`huron_vidal`.

## Models

**Residual Helmholtz models**

- Cubic EoS: SRK, PR76, PR78, RKPR, PSRK
  - Alpha functions: Soave, RKPR, Mathias-Copeman
  - Mixing rules: QMR, QMRTD, CMR, CMRTD, MHV, HV, HV-NRTL, sDDLC
- GERG-2008 (multifluid)
- PC-SAFT

**Excess Gibbs models**

- NRTL, NRTL-HV
- UNIQUAC
- UNIFAC VLE, UNIFAC Dortmund, PSRK-UNIFAC

## Properties and phase equilibria

- Volume, pressure, fugacity coefficients and their derivatives
- Residual properties (H, S, G, U, Cv, Cp); activity coefficients and excess
  properties from Ge models
- Flash calculations at PT and VT, and multiphase flash
- Saturation points: bubble, dew and liquid-liquid; pure-component saturation
- Stability analysis (tangent plane distance)
- Phase envelopes: PT (isopleths), Px (isotherms) and Tx (isobars), including
  three-phase and general multiphase envelopes
- Critical points and critical lines
- Global phase equilibrium diagrams of binary mixtures (`GPEC`, Python)

## Adding new models

New models only need their energy function; the rest of the library works on
top of it. Derivatives can come from three routes:

- **Analytical** derivatives, coded by hand.
- **Hyperdual numbers**, using the `hyperdual` module developed by
  [Philipp Rehner](https://github.com/prehner) and
  [Gernot Bauer](https://github.com/g-bauer). A complete PR76 implementation
  is in [`example/extra/adiff_pr76.f90`](example/extra/adiff_pr76.f90).
- **Tapenade**-generated derivatives; see
  [`tools/tapenade_diff/`](tools/tapenade_diff/).

The [new models guide](https://ipqa-research.github.io/yaeos/page/usage/newmodels/index.html)
walks through each route.

## Other EoS libraries

yaeos is one of several open-source equation of state libraries, each with its
own design choices:

- [Clapeyron](https://github.com/ClapeyronThermo/Clapeyron.jl) (Julia)
- [FeOs](https://github.com/feos-org/feos) (Rust, with Python bindings)
- [teqp](https://github.com/usnistgov/teqp) (C++, with Python bindings)
- [thermo](https://github.com/CalebBell/thermo) (Python)
- [thermopack](https://github.com/thermotools/thermopack) (Fortran, with Python bindings)
- [CoolProp](https://github.com/CoolProp/CoolProp) (C++, with Python bindings)

yaeos emphasizes readable, extensible modern Fortran so that researchers can
implement new models without dropping to lower-level code, while keeping
good performance and robust phase-diagram algorithms.

## Contributing

Contributions are welcome. Read [CONTRIBUTING.md](CONTRIBUTING.md) and the
[developer handbook](https://ipqa-research.github.io/yaeos/page/developer/index.html)
before opening a pull request. Bug reports and questions go to the
[issue tracker](https://github.com/ipqa-research/yaeos/issues).

If you use VS Code, predefined settings for developing and debugging are
available; from the repository root run:

```bash
git clone https://github.com/ipqa-research/vscode-fortran .vscode
```

## Citing, license and versioning

yaeos is maintained by the research group of Prof. Martín Cismondi-Duarte at
[IPQA (UNC-CONICET)](https://ipqa.unc.edu.ar/en/).

- **Citing:** if yaeos is useful in your work, please cite it using
  [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button
  exports BibTeX and APA).
- **License:** [Mozilla Public License 2.0](LICENSE).
- **Versioning:** [Semantic Versioning](https://semver.org/).
