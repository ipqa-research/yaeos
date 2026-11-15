---
title: C interface layer
---

[TOC]

# The `yaeos_c` interface module

`c_interface/yaeos_c.f90` is the flat, array-based façade over the
polymorphic Fortran core. It exists so that a caller which cannot speak
Fortran — today that is `f2py`, used by the Python package — can build a
thermodynamic model and evaluate properties by passing plain numbers and
flat arrays.

This page is for developers changing that layer. If you are *calling*
`yaeos` from Python, read [the Python interface page](python_interface.html).
If you are adding a model, read
[adding your own models](../usage/newmodels/index.html) first. A model used only
by native Fortran callers need not change this layer, but exposing a new model
to Python normally requires constructor/registration support here. Read the file itself before editing
it — the file is the contract, this page only describes it.

## What the layer contains

One module, `yaeos_c`, with three groups of routines:

* **Instantiation** — `pr76`, `pr78`, `srk`, `rkpr`, `psrk`, `pcsaft`,
  `multifluid_gerg2008`, `nrtl`, `uniquac`, `unifac_vle`, `unifac_psrk`,
  `unifac_dortmund`. Builds the model, stores it in the registry, returns an
  integer handle.
* **Configuration** — `set_qmr`, `set_qmrtd`, `set_mhv`, `set_hv`,
  `set_hvnrtl`, `set_sddlc`, `set_cmr`, `set_cmrtd`, `set_alpha_*`,
  `pcsaft_set_kij`. Attaches a mixing rule, an alpha function, or interaction
  parameters to an existing handle.
* **Evaluation** — `pressure`, `volume`, `lnphi_vt`,
  `helmholtz_residual_pt`, `flash`, `saturation_pressure`,
  `pt_mp_phase_envelope`. Computes one property or one equilibrium for a
  given handle.

The module is `private` by default and re-exports its surface through an
explicit `public ::` list near the top of the file. That list is not
decoration: `f2py` only generates wrappers for names that are public, so a
new routine that is not listed is unreachable from Python. Units are the
library defaults, stated in the module header: pressure in bar, volume in
litre, temperature in kelvin, amount in mole.

## `ISO_C_BINDING` kinds are not a C ABI

This is the most important thing to understand before editing this file.

The module header contains one import:

```fortran
use iso_c_binding, only: c_double, c_int, c_int64_t
```

These kinds correspond to the associated C types: `real(c_double)` for
double-precision reals, `integer(c_int)` for handles and
counts. It does **not** export a symbol name, does **not** define a C
calling convention, and does **not** make scalar arguments pass by value.

`bind(c)` is **absent** from this interface module. The hyperdual type in
`src/adiff/hyperdual.f90` does use `bind(c)`, but a type declaration does not
create C-callable wrapper procedures. No procedure in `yaeos_c` carries
`bind(c)`. The
consequences are concrete: compilers mangle Fortran names, so
the module procedure `pr76` is not thereby a C symbol named `pr76`. Ordinary
Fortran assumed-shape arrays can require compiler-specific descriptors, and
character arguments can involve hidden length arguments. C-compatible kinds
do not specify these calling conventions. The current repository does not
provide a portable public C header for these procedures. Do not assume that
declaring matching-looking C prototypes makes direct calls portable.

Some background material about this project describes the layer as a
`BIND(C)` interface. Treat that text as design intent; the shipped file is
the contract, and it is an `f2py`-oriented façade with C-compatible *kinds*,
not a stable C ABI. If you genuinely need a C entry point, do not retrofit
`bind(c)` onto `yaeos_c.f90` without designing and testing the ABI. Decide how
arrays, strings, optional values, handles, and errors cross the boundary.
For a simple explicit-shape ABI, pass extents explicitly and use interoperable
character buffers or lengths. Scalars may be passed by reference, or use
`value` to match a by-value C prototype. Modern Fortran also supports additional
interoperability mechanisms, including C descriptors and some optional
arguments; these require matching headers and compiler support, not guesses.
A small separate adapter with documented prototypes and C-side tests can
preserve the existing f2py layer while establishing a deliberate C API.

## Handles: the model registry

The core models are polymorphic objects (`class(ArModel)`,
`class(GeModel)`), and Fortran cannot expose polymorphic objects through a
flat interface. The layer therefore keeps a registry and hands out integer
handles.

```fortran
type :: ArModelContainer
   class(ArModel), allocatable :: model
end type ArModelContainer

type :: GeModelContainer
   class(GeModel), allocatable :: model
end type GeModelContainer

integer, parameter :: max_models  = 1000000 !! Maximum number of models
logical :: free_ar_model(max_models) = .true.
logical :: free_ge_model(max_models) = .true.
class(ArModelContainer), allocatable :: ar_models(:)
class(GeModelContainer), allocatable :: ge_models(:)
```

`ar_models` and `ge_models` are allocated lazily, on the first model
creation, with `max_models` entries each. The handle is a one-based index
into that table: `id = 1` is the first slot.

Instantiation follows a fixed recipe — build into a module-level singleton,
then move the singleton into the first free slot:

```fortran
subroutine pr76(tc, pc, w, id)
   use yaeos, only: PengRobinson76
   real(c_double), intent(in) :: tc(:)
   real(c_double), intent(in) :: pc(:)
   real(c_double), intent(in) :: w(:)
   integer(c_int), intent(out) :: id

   ar_model = PengRobinson76(tc, pc, w)
   call extend_ar_models_list(id)
end subroutine pr76
```

`extend_ar_models_list` scans `free_ar_model` from index 1, marks the first
free slot as taken, and `move_alloc`s the singleton into it. The `GeModel`
path is identical through `extend_ge_models_list`.

### Lifetime and release

`make_available_ar_models_list(id)` sets `free_ar_model(id) = .true.` and
nothing else. It does **not** deallocate the model; the object stays in the
container until the slot is reused by a later `move_alloc`. Consequently,
releasing a handle returns the *slot*, not the memory, and a released
handle is not invalidated — using it afterwards is undefined behaviour,
because the slot may already hold a different model. Treat handles as
opaque and never cache them beyond the object that owns them.

Because the registry and the singleton buffers are module state, the module
is not re-entrant. Do not create models concurrently from several threads.
Do not infer that evaluation is thread-safe either: inspect the called model,
solver, AD runtime, and any shared flags or caches. Read-only handle lookup
alone is insufficient to establish reentrancy of the entire calculation.

The exhaustion guard deserves a note:

```fortran
if (id == max_models) error stop 1
```

`id` is `intent(out)` and is only assigned inside the scan loop, so when
the table is full this comparison inspects a value the procedure never set.
Do not rely on this branch to produce a clean error. The check also rejects
the final slot after allocation. Avoid registry exhaustion and test any future
capacity/error-handling change explicitly.

## Argument conventions a caller must respect

**Order, names and returns.** Dummy arguments are positional. `f2py` also
exposes them as keyword arguments, lowercased (`--lower` is passed to
`f2py` in the Meson build), so the Fortran name `dPdn` becomes the Python
keyword `dpdn`. An `intent(out)` scalar becomes the value of the call: in
`pressure` the pressure `P` is `intent(out)`, so `f2py` returns it and the
remaining optional outputs become keyword arguments.

**`optional, intent(inout)` arguments are optional keywords.** Omit the
keyword and the generated binding can represent its absence; wrapper code
uses presence checks to decide which optional outputs to pass onward. Verify
this behavior in the generated binding when changing an argument. Shared
property kernels may still compute derivative intermediates even when no
derivative is returned to the user.

**Shapes.** Several dummy arguments take their extent from another
argument:

```fortran
real(c_double), optional, intent(in out) :: dPdV
real(c_double), optional, intent(in out) :: dPdT
real(c_double), optional, intent(in out) :: dPdn(size(n))
real(c_double), intent(out) :: x(size(z0))
real(c_double), intent(out) :: x_ls(max_points, np, size(z))
```

A direct caller must supply matching extents. `f2py` emits parser
diagnostics for some dependent dimensions; those diagnostics alone do not
establish what shape a generated wrapper accepts or returns. Check the
generated signature and exercise the relevant call when changing a shape.

Some outputs have a fixed size chosen by this file rather than by the
caller: `pure_saturation_line` fills arrays of 10000 elements, and
`find_self_intersections` fills arrays of 10.

**Character arguments are fixed-width and padded.** Several `root_type` and
`kind` arguments use `character(len=15)`. Root selectors use `"stable"`,
`"liquid"`, and `"vapor"`; saturation-point kinds instead use values such as
`"bubble"`, `"dew"`, and `"liquid-liquid"`. Follow the specific procedure's
contract. `convert_kind` maps the numeric root selectors used by multiphase
envelopes back to root strings.

**Integer kinds are mixed.** Most arguments use `integer(c_int)`, but
several use the plain default `integer`: `nc`, `ngs`, `g_ids` and `g_v` in
`psrk`, `ids` in `multifluid_gerg2008`, and `spec`, `max_iters` and
`max_points` in the critical-point routines. `f2py` maps both to a C `int`
in the generated configurations used by common builds, but do not assume
default integer always matches `c_int`. Check the generated type mapping when
changing compiler kinds or portability settings; normalizing kinds is a
separate compatibility change that needs tests.

## Column-major layout

Fortran arrays are column-major. For a two-dimensional array the element
at Fortran position `(i, j)` sits at zero-based C offset
`(i - 1) + (j - 1) * nrows`.

This matters whenever a matrix crosses the boundary — the NRTL `aij`,
`bij`, `cij` matrices, the mixing-rule `kij`/`lij` matrices, the group
matrices of the UNIFAC constructors, or a Hessian of mole derivatives. A
future explicit-buffer C adapter would need to document this layout. `f2py` is more
forgiving, but only partly, and the difference decides whether a wrong
answer is loud or silent:

* An `intent(in)` two-dimensional argument is converted to Fortran order
  by `f2py` when the caller's array is not already Fortran-contiguous.
  The logical values are preserved; only a copy is added.
* An `intent(inout)` two-dimensional argument is **not** converted. The
  caller must supply a Fortran-contiguous array, and `f2py` raises
  `ValueError: failed to initialize intent(inout) array -- input not
  fortran contiguous` otherwise. This is why every Hessian and
  derivative matrix in `python/yaeos/core.py` is allocated with
  `order="F"`.
* A dummy argument **without** an `intent` attribute is not equivalent to
  explicitly declaring `intent(inout)`. Fortran imposes no intent restriction,
  and f2py may infer or default to its own treatment. `set_qmr` (`kij`, `lij`)
  and `set_cmr` (`kijk`, `lijk`) include such arguments. Check the generated
  binding rather than assuming that writes propagate back to the caller;
  a temporary input copy need not provide write-back behavior.

The Python layer avoids the last case by always passing an
`np.array(matrix, order="F")` copy, as in
`python/yaeos/models/excess_gibbs/nrtl.py` and in
`python/yaeos/models/residual_helmholtz/cubic_eos/mixing_rules.py`.

## Optional outputs and NaN sentinels

Two different conventions express "no value", and they are not
interchangeable.

**NaN marks missing, invalid, or unfilled output, depending on the routine.** A helper function computes
`0.0 / 0.0` and is used to pre-fill output arrays before the points that
were actually computed are written in. A reader must therefore filter on
NaN rather than assume that a returned array is dense. `critical_line`
fills `as`, `Ts`, `Ps` and `Vs` with NaN and sets every CEP output to NaN
when the critical line has no critical end point; `binary_llv_from_cep`
fills its eight output arrays with NaN; the three multi-phase envelope
routines pre-fill `x_ls`, `ws`, `betas`, `Ps`, `Ts`, `alphas`, `Pcs` and
`Tcs`, so trailing points of a short envelope read as NaN rather than as
zeros; `pure_saturation_line` pre-fills its four 10000-element arrays and
writes only the computed prefix; `find_self_intersections` pre-fills its
outputs.

**Failure sentinels are routine-specific.** When a saturation
calculation returns a state whose `kind` is `"failed"`, the wrapper sets
the pressure or the temperature to NaN instead of returning a number. The
flash wrappers use a different signal: `beta = -1`, with `x` and `y` set
to the feed composition and `Vx` and `Vy` set to `1`. `flash_second_order`
additionally checks `trim(result%kind) == "failed"`, which plain `flash`
does not. Check for failure before using a result: a NaN propagates
through NumPy arithmetic without raising, and `beta == -1` is easy to miss
in a dictionary.

## Trace of one call: pressure

**Step 1 — the flat wrapper** in `c_interface/yaeos_c.f90` takes the
handle, the amounts, the volume and the temperature, returns the pressure
through an `intent(out)` argument, and offers the three derivatives as
`optional, intent(inout)` arguments:

```fortran
subroutine pressure(id, n, V, T, P, dPdV, dPdT, dPdn)
   integer(c_int), intent(in) :: id
   real(c_double), intent(in) :: n(:)
   real(c_double), intent(in) :: V
   real(c_double), intent(in) :: T
   real(c_double), intent(out) :: P
   real(c_double), optional, intent(in out) :: dPdV
   real(c_double), optional, intent(in out) :: dPdT
   real(c_double), optional, intent(in out) :: dPdn(size(n))

   call ar_models(id)%model%pressure(n, V, T, P, dPdV, dPdT, dPdn)
end subroutine pressure
```

**Step 2 — dynamic dispatch.** `ar_models(id)%model` is the polymorphic
object stored in the registry slot. The `select type` blocks elsewhere in
the file, for instance in the mixing-rule setters, show the same pattern:
the concrete type is resolved at run time.

**Step 3 — the model method** in
`src/models/residual_helmholtz/ar_models.f90` calls `residual_helmholtz`
for the energy and its derivatives and then adds the ideal-gas term
itself: `P = totn * R * T/V - ArV`, with `dPdV = -ArV2 - R*T*totn/V**2`,
`dPdT = -ArTV + totn*R/V` and `dPdn(:) = R*T/V - ArVn(:)`. The ideal-gas
term is added by the inherited `pressure` method, not by the residual
Helmholtz energy; do not add it again in a wrapper. Note also the kind
change: the model works in `real(pr)`, and `pr` is `real64`
(`src/constants.f90`). Verify that the chosen compiler's `pr` and `c_double`
are compatible when porting the build; do not rely only on similar kind names.

## Adding a routine to this layer

1. Write the subroutine in `c_interface/yaeos_c.f90`, inside the module,
   with `iso_c_binding` kinds on every dummy argument.
2. Add its name to the `public ::` list. Without this, `f2py` will not
   generate a wrapper and Python will not see it.
3. Use `intent(in)` for inputs and `intent(out)` for results. Mark
   derivative outputs `optional, intent(inout)` so callers can request
   only what they need.
4. Take the model from the registry: `ar_models(id)%model` or
   `ge_models(id)%model`. Never cache a model in a local variable across
   calls.
5. Decide the failure convention explicitly: return NaN for the property
   or the documented sentinel, and pre-fill array outputs with NaN when
   the number of meaningful entries can be smaller than the array extent.
   Say both in a comment.
6. Add the Python wrapper that allocates the derivative buffers and shapes
   the result. See [the Python interface page](python_interface.html).
7. Add a test. The Python-side test suite lives in `python/tests/`, and
   `python/tests/models/residual_helmholtz/test_ar_same_as_fortran.py`
   compares the Python result against a Fortran program it generates and
   runs. This checks the interface against the native implementation; add
   independent reference or mathematical checks for the property itself.
8. Keep the doc comment style of the file. The `!!` comments follow the
   FORD convention, but `ford.md` sets `src_dir: src`, so this directory is
   **not** part of the generated API documentation. The comments are still
   worth writing for IDE users and maintainers of these routines.

## Debugging boundaries

| Symptom | Look first at | Then at |
|:--|:--|:--|
| `AttributeError` on `yaeos_c` | The `public ::` list in `yaeos_c.f90` | Whether the extension was rebuilt |
| Wrong numbers for a matrix argument | The `intent` of the dummy argument | Whether the routine writes into the matrix |
| Wrong numbers for a vector argument | Component ordering against the model parameters | Units: bar, litre, kelvin, mole |
| A derivative is garbage but the value is fine | The optional buffer was not allocated F-ordered | Whether the model supports that derivative |
| `NaN` where a number was expected | The failure convention of that routine | The underlying equilibrium solver |
| `beta == -1` from a flash | The feed composition and the specification | Whether an initial guess was supplied |
| A handle refers to another model | A handle was released and reused | The lifetime rules above |
| The process aborts during model creation | The registry exhaustion guard | The number of live models |

## Known sharp edges

* No `bind(c)`, so this is not a C ABI and there is no header.
* `make_available_*_models_list` releases the slot, not the memory.
* The registry exhaustion guard reads an unassigned `intent(out)` argument.
* `c_int64_t` is imported and never used, and integer kinds are mixed
  between `integer(c_int)` and the default `integer`.
* `character(len=15)` is padded and must be trimmed by the caller.
* `flash_grid` prints progress to standard output when its `parallel`
  argument is false.
* The registry tables hold a million entries and are allocated on first use.
* The `!!` doc comments here are not rendered by FORD.

## Related pages

* [Developer handbook](index.html) — the rest of the developer guide.
* [Architecture](architecture.html) — how the Fortran layers fit together.
* [Fortran internals](fortran.html) — the core library this layer wraps.
* [Numerics](numerics.html) — the solvers behind the equilibrium routines.
* [Development workflow](workflow.html) — how to build and test a change.
* [Testing](testing.html) — the test suites and how to extend them.
* [Python interface layer](python_interface.html) — the wrapper above this file.
* [Adding your own models](../usage/newmodels/index.html) — the user-facing guide for extending the library.
