---
project: yaeos
summary: <div class="yaeos-hero">
         <img class="yaeos-hero-logo" src="media/logo.png" alt="yaeos logo">
         <div>
         <p class="yaeos-hero-tagline">Thermodynamic properties and phase
         equilibria with equations of state, in modern Fortran and Python.</p>
         <p class="yaeos-hero-sub">Robust algorithms for phase diagrams and
         an easy framework for adding new models, with analytical derivatives
         or automatic differentiation.</p>
         <div class="yaeos-hero-actions">
         <a class="btn btn-primary" href="page/usage/index.html">Get started</a>
         <a class="btn btn-outline-primary" href="module/yaeos.html">Fortran API</a>
         <a class="btn btn-outline-primary" href="page/python-api/index.html">Python API</a>
         <a class="btn btn-outline-secondary" href="https://github.com/ipqa-research/yaeos">GitHub</a>
         <a class="btn btn-outline-secondary" href="https://colab.research.google.com/github/ipqa-research/yaeos/blob/main/notebooks/run_online.ipynb">Try it in Colab</a>
         </div>
         </div>
         </div>
project_github: https://github.com/ipqa-research/yaeos
author: Federico E. Benelli
author_description: PhD student with focus on reservoir PVT simulation.
email: federico.benelli@unc.edu.ar
github: https://github.com/fedebenelli
src_dir: src
exclude_dir: test doc src/adiff/autodiff_api/tapenade/c src/adiff/autodiff_api/tapenade/PUSHPOPGeneralLib
output_dir: doc/ford_site
preprocessor: gfortran -E
display: public
         protected
         private
source: false
proc_internals: true
sort: permission-alpha
docmark_alt: !>
docmark: !
predocmark_alt: *
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
               markdown.extensions.md_in_html
graph: true
license: MPL
page_dir: doc/page
media_dir: doc/media
css: doc/css/yaeos.css
favicon: doc/media/logo.png
---

<div class="yaeos-features">
<div class="yaeos-feature">
<h3>Phase diagrams that close</h3>
<p>Continuation methods trace PT, Px and Tx envelopes through critical points,
including three-phase and multiphase boundaries, critical lines and binary
global phase diagrams.</p>
</div>
<div class="yaeos-feature">
<h3>Models you can extend</h3>
<p>Write the residual Helmholtz or excess Gibbs energy and get fugacities,
residual properties and every equilibrium routine for free. Derivatives can
be analytical, hyperdual or generated with Tapenade.</p>
</div>
<div class="yaeos-feature">
<h3>Fortran speed, Python reach</h3>
<p>A readable modern Fortran core built with <code>fpm</code>, wrapped as a
<code>pip install yaeos</code> package for notebooks, fitting and plotting.</p>
</div>
</div>

## Quick start

The same calculation in both languages: build a Peng-Robinson model for
methane + n-butane, evaluate a pressure, a bubble point and the full PT phase
envelope. Units are K, bar, L and mol throughout.

<ul class="nav nav-tabs yaeos-tabs" role="tablist">
<li class="nav-item" role="presentation"><button class="nav-link active" id="tab-fortran" data-bs-toggle="tab" data-bs-target="#pane-fortran" type="button" role="tab" aria-controls="pane-fortran" aria-selected="true">Fortran</button></li>
<li class="nav-item" role="presentation"><button class="nav-link" id="tab-python" data-bs-toggle="tab" data-bs-target="#pane-python" type="button" role="tab" aria-controls="pane-python" aria-selected="false">Python</button></li>
</ul>
<div class="tab-content yaeos-tab-content" markdown="1">
<div class="tab-pane fade show active" id="pane-fortran" role="tabpanel" aria-labelledby="tab-fortran" markdown="1">

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

</div>
<div class="tab-pane fade" id="pane-python" role="tabpanel" aria-labelledby="tab-python" markdown="1">

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

</div>
</div>

Both print a pressure of 15.91 bar, a bubble pressure of 57.37 bar and a
critical point near 401.5 K and 68.1 bar. The
[guide](page/usage/index.html) walks through each step and explains how to
read the results.

## Install

<div class="yaeos-install" markdown="1">
<div markdown="1">

**Fortran** — add yaeos as an [fpm](https://fpm.fortran-lang.org)
dependency (BLAS and LAPACK are required, e.g.
`sudo apt install libblas-dev liblapack-dev`):

```toml
[dependencies]
yaeos = {git="https://github.com/ipqa-research/yaeos"}
```

</div>
<div markdown="1">

**Python** — binary wheels are published on PyPI:

```bash
pip install yaeos
```

</div>
</div>

## Models

| Residual Helmholtz (\(A^r\)) | Excess Gibbs (\(G^E\)) |
|------------------------------|------------------------|
| [Cubic EoS](page/usage/eos/cubics/index.html): SRK, PR76, PR78, RKPR, PSRK | [NRTL](page/usage/excessmodels/nrtl.html) and [NRTL-HV](page/usage/excessmodels/nrtl-hv.html) |
| [Alpha functions](page/usage/eos/cubics/alpha.html): Soave, RKPR, Mathias-Copeman | [UNIQUAC](page/usage/excessmodels/uniquac.html) |
| [Mixing rules](page/usage/eos/cubics/mixing.html): QMR, QMRTD, CMR, CMRTD, MHV, HV, HV-NRTL, sDDLC | [UNIFAC VLE](page/usage/excessmodels/unifac_models/unifaclv.html), [Dortmund](page/usage/excessmodels/unifac_models/dortmund.html) and [PSRK](page/usage/excessmodels/unifac_models/psrk.html) |
| [GERG-2008](page/usage/eos/gerg2008.html) multifluid | |
| [PC-SAFT](page/usage/eos/pc_saft.html) | |

Missing one? [Add your own model](page/usage/newmodels/index.html) by
writing its energy function.

## Explore the documentation

<div class="yaeos-cards">
<a class="yaeos-card" href="page/usage/fundamentals.html"><h3>Fundamentals</h3><p>Equations of state, fugacity and phase equilibrium for programmers who are new to thermodynamics.</p></a>
<a class="yaeos-card" href="page/usage/eos/index.html"><h3>Equations of state</h3><p>Volumes, pressures, fugacities and residual properties from \(A^r\) models.</p></a>
<a class="yaeos-card" href="page/usage/excessmodels/index.html"><h3>Excess Gibbs models</h3><p>Activity coefficients and excess properties from \(G^E\) models.</p></a>
<a class="yaeos-card" href="page/usage/phase_equilibrium/index.html"><h3>Phase equilibrium</h3><p>Flashes, saturation points, stability, envelopes and critical points.</p></a>
<a class="yaeos-card" href="page/usage/newmodels/index.html"><h3>New models</h3><p>Implement a model with analytical derivatives or automatic differentiation.</p></a>
<a class="yaeos-card" href="page/developer/index.html"><h3>Developer handbook</h3><p>Architecture, numerics, the Python bridge, testing and workflow.</p></a>
<a class="yaeos-card" href="page/python-api/index.html"><h3>Python API</h3><p>Reference for the <code>yaeos</code> Python package, with tutorials.</p></a>
<a class="yaeos-card" href="lists/procedures.html"><h3>Fortran API</h3><p>Every module, type and procedure, generated from the source.</p></a>
<a class="yaeos-card" href="page/contributing/index.html"><h3>Contributing</h3><p>How to report issues, propose changes and follow the style guide.</p></a>
</div>

## Citing and community

yaeos is developed by the research group of Prof. Martín Cismondi-Duarte at
[IPQA (UNC-CONICET)](https://ipqa.unc.edu.ar/en/) and released under the
[Mozilla Public License 2.0](https://github.com/ipqa-research/yaeos/blob/main/LICENSE).
If yaeos is useful in your work, please cite it using the repository's
[`CITATION.cff`](https://github.com/ipqa-research/yaeos/blob/main/CITATION.cff)
(GitHub's "Cite this repository" button exports BibTeX and APA). Questions,
bug reports and ideas are welcome on the
[issue tracker](https://github.com/ipqa-research/yaeos/issues).
