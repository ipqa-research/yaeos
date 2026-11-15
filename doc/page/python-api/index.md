---
title: Python-API
---

[TOC]

# Python API

The Python package exposes the same models and phase-equilibrium routines as
the Fortran library, so you can build a model, evaluate properties, and
solve flashes from a Python script or notebook.

Start with the [Python tutorial](https://ipqa-research.github.io/yaeos/page/python-api/tutorial/tutorial_tree.html)
and use the API reference for exact argument and return-value details.

## Choose the right starting point

| You want to know | Start with |
|:--|:--|
| The pressure or density of a fluid at given conditions | An [equation of state](../usage/eos/index.html) |
| How nonideal a liquid mixture is (activity coefficients) | An [excess Gibbs model](../usage/excessmodels/index.html) |
| Whether a mixture separates into phases, and how much of each | A [phase-equilibrium calculation](../usage/phase_equilibrium/index.html) |
| The equations, derivatives, and references behind a model | The [Fortran user documentation](../index.html) |
| How the Python package is built and tested | The [developer handbook](../developer/index.html) |

## Conventions that prevent silent mistakes

- The property methods use **kelvin**, **bar**, **litres**, and **moles**;
  energies are reported in bar L.
- `moles` is an amount vector, while `z` in equilibrium calculations is an
  overall mole-fraction vector that must sum to one.
- Keep the component order identical in every parameter array, composition
  vector, and interaction matrix. The order used when the model was built is
  the order the results come back in.
- Results named `residual` or `excess` are departures from a reference state,
  not complete enthalpies or Gibbs energies.

## About this page

This page is a placeholder for a FORD-only documentation build. In CI, the
Sphinx-generated Python documentation is copied over the whole
`page/python-api/` directory, so this page and its `index.html` are replaced
by the Python tutorial and API reference listed above. If you are reading
this page in a local FORD build, the Python reference has not been copied
into that location yet.

The sources of the Python documentation live in `python/docs/source/`. Build
them with:

```bash
cd python/docs
sphinx-build -b html source build
```

The result is a self-contained HTML tree whose entry point is
`build/index.html`; copy it to `doc/ford_site/page/python-api/` to reproduce
what CI publishes.
