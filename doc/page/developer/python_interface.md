---
title: Python interface layer
---

[TOC]

# The `yaeos` Python package

The Python package is a thin, typed layer over the Fortran library. It
gives a Python user an object-oriented API — models with methods, results
as dictionaries, derivatives as named entries — instead of a flat list of
functions that take and return bare arrays.

This page is for developers changing that layer. For *using* the package,
read the [Python API guide](../python-api/index.html). For the layer that
sits underneath, read [the C interface page](c_interface.html). For
adding a thermodynamic model, read
[adding your own models](../usage/newmodels/index.html) first. Native Fortran
extension does not automatically publish a Python class: a new model exposed
to Python normally needs wrapper construction, a handle-owning class, and
package/build exports, even though shared property methods can be inherited.

## Package map

Everything lives under `python/`.

* `python/pyproject.toml` — package metadata and the `meson-python` backend.
* `python/meson.build` — the real build: fpm, f2py, and the extension module.
* `python/yaeos/__init__.py` — public exports; `__version__` from the metadata.
* `python/yaeos/core.py` — the `ArModel` and `GeModel` abstract classes and
  every property method.
* `python/yaeos/constants.py` — `R` and the `root_kinds` mapping.
* `python/yaeos/envelopes.py` — `PTEnvelope`, `PXEnvelope`, `TXEnvelope`.
* `python/yaeos/gpec.py` — the generalised phase-envelope calculator.
* `python/yaeos/models/` — concrete models: `excess_gibbs/`,
  `residual_helmholtz/`, `groups.py`.
* `python/yaeos/fitting/` — parameter-fitting helpers and their solvers.
* `python/yaeos/tools/writers.py` — writers that render model parameters as
  Fortran source.
* `python/yaeos/lib/__init__.py` — the boundary to the compiled extension.
* `python/tests/` — the Python test suite; `python/docs/source/` — Sphinx.

`python/yaeos/core.py` is large because it holds every thermodynamic
property method for both abstract classes; the concrete model classes
mostly hold parameters and call `yaeos_c`.

## The build pipeline: fpm, then f2py, then Meson

The full build instructions, the toolchain list and the working-directory
pitfalls are in the [developer workflow](workflow.html). What matters for
the Python layer is the order of the three stages, because it decides what
has to be rebuilt after a change.

**Stage 1 — the Fortran core.** `python/meson.build` shells out to
`fpm install --profile release --prefix <build>/fpm_build`, which produces
the static library `libyaeos.a` and the module files. The `include`
subdirectory matters: Meson adds `fpm_build/include` to the include path
so that `use yaeos` resolves when the extension is compiled.

**Stage 2 — the f2py wrapper.** A `custom_target` runs f2py over the flat
interface module:

```
python -m numpy.f2py @INPUT@ -m yaeos_python --lower
```

`@INPUT@` is `c_interface/yaeos_c.f90`, resolved as
`meson.global_source_root() / '..' / 'c_interface' / 'yaeos_c.f90'`. Two
files come out: `yaeos_pythonmodule.c` and
`yaeos_python-f2pywrappers2.f90`. The `--lower` flag lowercases the
Fortran names, which is why the Python side calls `yaeos_c.cv_residual_vt`
and `yaeos_c.set_alpha_rkpr` for Fortran routines named `Cv_residual_vt`
and `set_alpha_RKPR`.

f2py prints a long list of messages like

```
analyzevars: could not parse dimension for variable 'size(n)'
```

These messages indicate that f2py could not determine the dimension while
parsing the source. They do not, by themselves, establish whether the
generated wrapper handles that argument correctly. Inspect the generated
signature and test the actual accepted and returned shapes when changing
dependent-dimension arguments.

**Stage 3 — the extension module.** Meson compiles `yaeos_c.f90`
together with the generated files into a single extension module named
`yaeos_python` and installs it under the `yaeos/lib` subdirectory of the
package. The pure-Python modules are installed separately with
`py.install_sources`.

### Compiler requirements

`python/meson.build` enforces them, and the gate is worth knowing before
you install the package:

* Meson `>= 0.64.0`, and a Fortran compiler added with
  `add_languages('fortran', native: false)` so that a cross build is
  possible.
* On every platform except Windows, the Fortran compiler **must be version
  14.1 or newer**; otherwise the build stops with
  `Fortran compiler version 14.1 or newer is required!`. The comparison is
  `>=14.1`, so a GCC 14.0 toolchain fails. The flags added for that path
  — `-ftrampoline-impl=heap`, `-ffree-line-length-none`,
  `-fno-strict-aliasing`, `-fno-common`,
  `-fno-inline-small-functions` — are GNU-style, not portable.
* Windows skips the version check and links statically instead.
* `quadmath` is required everywhere except Windows, where it is linked
  explicitly, and OpenMP is linked with `-lgomp`.
* The build requires `fpm` and `ninja` on non-Windows platforms, declared
  as build requirements in `python/pyproject.toml`.

The Fortran core has no such requirement: `fpm` builds the library and its
tests with whatever compiler you give it. The 14.1 gate lives only in the
Python package build.

One inconsistency to know about: `python/meson.build` declares
`version : '3.1.0'`, while `python/pyproject.toml` and `fpm.toml` both
say `4.5.4`. The version the user sees comes from the installed
distribution metadata, so the Meson project version is stale.

### What needs a rebuild

| Change | Rebuild needed? |
|:--|:--|
| A pure-Python file under `python/yaeos/` | No |
| `c_interface/yaeos_c.f90` | Yes — the extension is compiled from it |
| A file under `src/` | Yes — `fpm install` re-runs |
| `python/meson.build` or `python/pyproject.toml` | Yes — Meson reconfigures |

## The `yaeos.lib` boundary

`python/yaeos/lib/__init__.py` is three lines long and defines the whole
boundary:

```python
from yaeos.lib.yaeos_python import yaeos_c

__all__ = ["yaeos_c"]
```

The compiled extension is `yaeos_python`; the Fortran module `yaeos_c`
inside it appears as the submodule `yaeos_c`. Every routine of the flat
interface is therefore reachable as
`yaeos_c.<fortran_name_lower_cased>`.

Application code should use the public classes exported by
`python/yaeos/__init__.py` — `PengRobinson76`, `NRTL`, `QMR` and so on.
The `yaeos_c` namespace is an implementation detail with no stability
guarantee. And nothing above `python/yaeos/lib/` may hand NumPy arrays to
`yaeos_c` with the wrong memory order or the wrong shape: that is where a
silent wrong answer enters the system.

## Trace: `PengRobinson76` to a pressure

The whole path of one number, from a Python constructor to a Fortran
residual Helmholtz evaluation and back.

**Step 1 — the constructor** in
`python/yaeos/models/residual_helmholtz/cubic_eos/cubic_eos.py` stores the
parameters and asks the flat layer for a handle:

```python
class PengRobinson76(CubicEoS):
    def __init__(self, critical_temperatures, critical_pressures,
                 acentric_factors, mixrule: CubicMixRule = None) -> None:
        super(PengRobinson76, self).__init__(
            critical_temperatures, critical_pressures, acentric_factors
        )
        self.id = yaeos_c.pr76(self.tc, self.pc, self.w)
        self.mixrule = mixrule
        if mixrule:
            mixrule.set_mixrule(self.id)
```

`self.id` is the handle used by this Python object's calls into the Fortran
registry. The model itself remains in module-level Fortran registry state.

**Step 2 — the handle registry.** `yaeos_c.pr76` calls the Fortran
`pr76`, which builds the concrete model, stores it in the first free slot
of the `ArModel` registry, and returns the one-based index of that slot.
The registry is described in [the C interface page](c_interface.html); the
point to remember here is that the handle is an opaque integer owned by
the Fortran module.

**Step 3 — the property wrapper** in `python/yaeos/core.py` allocates the
output buffers and calls the flat layer:

```python
def pressure(self, moles, volume, temperature, dv=False, dt=False, dn=False):
    nc = len(moles)

    dv = np.array(0, dtype=np.float64) if dv else None
    dt = np.array(0, dtype=np.float64) if dt else None
    dn = np.empty(nc, order="F") if dn else None

    res = yaeos_c.pressure(
        self.id, moles, volume, temperature, dpdv=dv, dpdt=dt, dpdn=dn
    )
```

Note the pattern: a boolean flag becomes an output buffer, and `None`
means "do not ask for this derivative". The buffers are allocated with
`order="F"` because the derivative matrices are Fortran-ordered on the
other side.

**Step 4 — the f2py marshalling.** `yaeos_c.pressure` reaches the
generated C wrapper, which converts the Python objects to Fortran
arguments, calls the Fortran routine, and converts the results back. The
rules that matter, each verified against a small `f2py` probe built from
this repository's interface file:

* An `intent(out)` scalar becomes the value of the call. The pressure `P`
  is `intent(out)`, so the call returns a float. Several `intent(out)`
  arguments come back as a tuple in declaration order: `flash` has seven
  of them, and the Python wrapper unpacks them as
  `x, y, pressure, temperature, volume_x, volume_y, beta`.
* An `optional, intent(inout)` argument becomes an optional keyword.
  Omitting it, or passing `None`, means the argument is not passed at all
  and the Fortran side sees `.not. present(...)`. A rank-0 argument wants
  a zero-dimensional NumPy array; `core.py` uses
  `np.array(0, dtype=np.float64)` and reads the result back with
  `float(...)`.
* A one-dimensional array may be in any order. A two-dimensional array
  that Fortran *writes into* must be Fortran-contiguous: passing a
  C-ordered one raises `ValueError: failed to initialize intent(inout)
  array -- input not fortran contiguous`. A two-dimensional array that
  Fortran only *reads* is converted by `f2py`, so its logical values are
  always right.
* A `character(len=15)` argument such as `root_type` receives a Python
  string blank-padded to fifteen characters.

**Step 5 — the Fortran model.** The wrapper calls
`ar_models(id)%model%pressure(...)`, which dispatches dynamically to the
concrete model, which evaluates the residual Helmholtz energy and adds the
ideal-gas term. This is where the physics happens; everything above it is
plumbing.

**Step 6 — the result.** Back in Python the wrapper reshapes the outcome
into a consistent shape: a bare value when no derivative was requested,
and a `(value, derivatives)` pair otherwise, with lowercase derivative
keys.

**Step 7 — the lifetime.** `ArModel.__del__` releases the Fortran slot
with `yaeos_c.make_available_ar_models_list(self.id)`; `GeModel.__del__`
does the same with `make_available_ge_models_list`. Releasing a slot
returns it to the registry; the model object itself is only released when
the slot is reused. `__del__` is best-effort, so a session that creates and
drops many models in a loop should not assume the memory comes back
immediately.

## Shapes, memory order and dtypes

| Fortran dummy argument | Python buffer | Rule |
|:--|:--|:--|
| `n(:)` amounts | 1-D `float64` | Order irrelevant |
| `dPdn(size(n))` | 1-D `float64` | Order irrelevant; `order="F"` is used for consistency |
| `Arn2(size(n), size(n))` | 2-D `float64` | **Must** be Fortran-contiguous |
| `aij(:,:)` and friends | 2-D `float64` | Read-only; `f2py` converts |
| `kij(:,:)` in `set_qmr` | 2-D `float64` | The Fortran dummy has no `intent`, so `f2py` may not convert it |
| `root_type` | `str` | One of `"stable"`, `"liquid"`, `"vapor"` |

The last two rows are the sharpest edges in this layer. In
`c_interface/yaeos_c.f90`, `set_qmr` declares `kij` and `lij` and
`set_cmr` declares `kijk` and `lijk` **without an `intent` attribute**.
A dummy argument without `intent` carries no intent restriction (which is
not the same as declaring `intent(inout)`), and `f2py` does not enforce the
memory order for it: in a probe built from
this interface file, a value written by Fortran into `kij(1, 2)` of a
C-ordered array was not visible at `kij[0, 1]` in NumPy. The Python
classes avoid the problem by always passing an `order="F"` copy:

```python
class QMR(CubicMixRule):
    def __init__(self, kij, lij) -> None:
        self.kij = np.array(kij, order="F")
        self.lij = np.array(lij, order="F")

    def set_mixrule(self, ar_model_id: int) -> None:
        yaeos_c.set_qmr(ar_model_id, self.kij, self.lij)
```

Follow that convention in any new model class: a C-ordered matrix handed
to one of those routines may be misinterpreted without an error being
raised. Note also that `NRTL.__init__` stores the Fortran-ordered copies
in `self.a`, `self.b` and `self.c` but passes the original arguments to
`yaeos_c.nrtl` — harmless there because the `nrtl` dummies are
`intent(in)`, so `f2py` converts them, but not the pattern to copy.

## Return conventions

| Kind of result | Shape |
|:--|:--|
| Scalar property | `float` |
| Property with derivatives | `(value, {"dt": ..., "dv": ..., "dn": ...})` |
| Flash | `dict` with `x`, `y`, `Vx`, `Vy`, `P`, `T`, `beta` |
| Phase envelope | `PTEnvelope`, `PXEnvelope` or `TXEnvelope` object |
| Stability analysis | `({"w": ..., "tm": ...}, {"tm": [...], "w": [...]})` |

Two failure conventions cross the boundary and must be checked. `beta == -1`
means the flash did not converge: the compositions are then the feed
composition and the volumes are placeholders. `NaN` means a value does not
exist — a saturation pressure that failed, or an envelope point that was
never computed. NumPy propagates `NaN` silently, so filter on it
explicitly. The envelope classes in `python/yaeos/envelopes.py` already
drop the `NaN` rows; a raw call to `yaeos_c` does not.

## The package mirrors the Fortran API

`python/yaeos/core.py` defines `GeModel` and `ArModel` as abstract
classes that mirror the Fortran `class(GeModel)` and `class(ArModel)`.
Each concrete model inherits one of them and implements `size()`, plus
`_model_params_as_str()` and `_model_params_declaration_as_str()`, which
render Fortran source that constructs the model. Those last two are not
decoration: the equivalence tests use them to generate a Fortran program,
run it with `fpm`, and compare its output with the Python result, so a new
model class that omits them cannot take part in that check.

## Adding a feature: wrapper, Python, tests

Work through these in order; each step depends on the previous one.

1. **Add the flat routine** in `c_interface/yaeos_c.f90`, and add its name
   to the module's `public ::` list so that `f2py` generates a wrapper.
2. **Add the Python wrapper** in `python/yaeos/core.py`, on the abstract
   class that owns the property. Allocate the output buffers, pass the
   handle, and reshape the result into the house style: a bare value, or a
   `(value, derivatives)` pair with lowercase derivative keys.
3. **Expose it publicly** if it belongs to a concrete model: add the
   method to the class in `python/yaeos/models/`, and export the class
   from `python/yaeos/__init__.py` if it is new.
4. **Add a test** under `python/tests/`. For a property, extend
   `python/tests/models/residual_helmholtz/test_ar_same_as_fortran.py` so
   that the Fortran and the Python values are compared. For a model,
   follow the neighbouring tests in the same directory. The
   [testing guide](testing.html) covers what a good regression looks like.
5. **Update the Sphinx sources** under `python/docs/source/` if the new
   method is part of the public API.
6. **Keep the styles**: `black` with a line length of 79, `flake8` with
   `N803`, `W503` and `E203` ignored, and `pydocstyle` with the NumPy
   convention.

## Testing

The Python suite lives in `python/tests/` and mirrors the source tree:
`models/residual_helmholtz/` for the equations of state,
`models/gibbs_excess/` for the activity models, `models/algorithm/` for
the stability analysis, `models/equilibria/` for flash calculations, and
`fitting/` for the parameter fitting.

`python/tox.ini` defines the environments. The default one runs, in order,
`fpm run ge_python_test_values`, `pip install .`, and
`pytest tests/ --import-mode=importlib`. The first command runs the
Fortran helper in `app/ge_python_test_values.f90`, which writes
`ge_test_vals.txt` in the working directory; that file is gitignored. No
test in the current tree reads it, so treat it as a generated fixture
rather than a required step.

The strongest test in the suite is the Fortran–Python equivalence check.
It builds a Fortran program from the model's own parameter source, writes
it into the repository's `app/` directory, runs it with `fpm run`, and
compares the printed values with the Python results:

```python
fortran_code += "program main\n"
fortran_code += "use yaeos\n\n"
fortran_code += "implicit none\n\n"
fortran_code += ar_model._model_params_declaration_as_str()
fortran_code += ar_model._model_params_as_str()
fortran_code += thermoprops_code
fortran_code += "end program main"
```

Note that this test writes `app/test_ar_program.f90` into the working
tree. Do not commit it. Other environments: `tox -e style` runs
`flake8`, `tox -e docstyle` runs `pydocstyle`, `tox -e docs` builds the
Sphinx pages, and `tox -e coverage` runs the suite with `coverage` and
fails below 50 %.

## Debugging boundaries

| Symptom | Look first at | Then at |
|:--|:--|:--|
| `ImportError` on `yaeos.lib.yaeos_python` | The extension build in `python/meson.build` | The `fpm install` step and the Fortran compiler version |
| `AttributeError` on a `yaeos_c` function | The `public ::` list in `c_interface/yaeos_c.f90` | Whether the extension was rebuilt |
| `ValueError: ... not fortran contiguous` | The `order="F"` on the allocated buffer | The Fortran dummy's `intent` |
| A matrix argument gives wrong numbers | Whether the routine *writes* into the matrix | Whether the dummy declares an `intent` |
| A derivative is wrong but the value is right | The buffer allocated for that derivative | Whether the model supports that derivative at all |
| A returned array has `NaN` in it | The failure convention of that routine | The equilibrium solver |
| The Python and Fortran results disagree | The `_model_params_as_str` output | The component ordering of the model |
| The build stops at the Fortran version check | The compiler version in `python/meson.build` | The flags, which are GNU-specific |

## Known sharp edges

* The `yaeos_c` namespace is generated from the Fortran module and has no
  stability guarantee. Use the public classes.
* The extension module is built from `c_interface/yaeos_c.f90`, so a
  change there needs a rebuild, not just a reinstall.
* The Meson project version (`3.1.0`) disagrees with the package version
  (`4.5.4`).
* `set_qmr` and `set_cmr` declare their matrix arguments without an
  `intent`, which makes them sensitive to the caller's memory order.
* `NRTL.__init__` passes the original arrays rather than the
  Fortran-ordered copies it stores.
* `__del__` releases the Fortran slot, not the model memory.
* The equivalence tests write generated Fortran programs into `app/`.

## Related pages

* [Developer handbook](index.html) — the rest of the developer guide.
* [Architecture](architecture.html) — how the Fortran layers fit together.
* [Developer workflow](workflow.html) — the build, the toolchain, and the review loop.
* [Fortran internals](fortran.html) — the core library this layer wraps.
* [Numerics](numerics.html) — the solvers behind the equilibrium routines.
* [Testing](testing.html) — what a good regression test looks like.
* [End-to-end change](walkthrough.html) — planning a feature across all the layers.
* [C interface layer](c_interface.html) — the flat layer underneath this one.
* [Python API guide](../python-api/index.html) — the user-facing guide.
* [Adding your own models](../usage/newmodels/index.html) — the user-facing guide for extending the library.
