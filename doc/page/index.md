---
title: User documentation
ordered_subpage: usage
                 python-api
                 developer
                 contributing
---

# Learn to use yaeos

<p align="center">
  <img src="../media/logo.svg" width="500" alt="yaeos logo">
</p>

`yaeos` is a Fortran library for calculating fluid properties and phase
equilibria. It includes equations of state (EOS) for gases and liquids, and
activity-coefficient models for non-ideal liquid mixtures.

You can use it to answer questions such as “What pressure does this mixture
have?”, “When does the first liquid drop appear?”, or “How much liquid and
vapor coexist at this temperature and pressure?” These are different
calculations, and this guide helps you choose the right one.

## Start with a working example

<div class="quick-links">
  <a href="usage/index.html">Your first calculation</a>
  <a href="usage/fundamentals.html">Thermodynamics essentials</a>
  <a href="usage/eos/cubics/index.html">Choose an EOS</a>
  <a href="usage/excessmodels/index.html">Liquid activity models</a>
  <a href="usage/phase_equilibrium/index.html">Phase equilibrium</a>
  <a href="python-api/index.html">Python users</a>
  <a href="developer/index.html">Become a developer</a>
</div>

The [getting-started tutorial](usage/index.html) builds a small program,
predicts a mixture's pressure, solves for volume, and explains how to check
the result. If terms such as *mole fraction*, *fugacity*, or *volume root* are
unfamiliar, keep [thermodynamics essentials](usage/fundamentals.html) alongside
the tutorial. You do not need to derive an EOS before running it.

## A suggested learning path

| Step | Read | What you will be able to do |
|------|------|-----------------------------|
| 1 | [Getting started](usage/index.html) | Install the library and run a complete Fortran calculation |
| 2 | [Thermodynamics essentials](usage/fundamentals.html) | Distinguish amounts, compositions, phases, units, and energy references |
| 3 | [EOS guides](usage/eos/cubics/index.html) or [activity models](usage/excessmodels/index.html) | Choose a model family and assemble its input data |
| 4 | [Property reference](usage/eos/index.html) | Evaluate properties and understand optional derivatives |
| 5 | [Phase equilibrium](usage/phase_equilibrium/index.html) | Select saturation, flash, stability, or envelope calculations |
| 6 | [Adding models](usage/newmodels/index.html) | Extend the library and check a new model's derivatives |
| 7 | [Developer handbook](developer/index.html) | Navigate the codebase and contribute across Fortran, C-facing wrappers, and Python |

## Choose a guide by the question you have

* **I know the state and want a property.** Start with the
  [property guide](usage/eos/index.html). EOS models share the `ArModel`
  interface, so property calls stay similar when you switch between models.
* **I have a non-ideal liquid mixture.** Start with
  [excess-Gibbs models](usage/excessmodels/index.html). Their `GeModel`
  interface supplies activity coefficients and excess properties, not a
  general pressure–volume relationship.
* **I want to know which phases exist.** Start with
  [phase equilibrium](usage/phase_equilibrium/index.html). A single-phase
  property call is not a phase-split calculation.
* **I use Python rather than Fortran.** Follow the
  [Python tutorials and API](python-api/index.html). The physical model and
  units are shared, but Python method names and result containers can differ.
* **I am extending or contributing to the project.** Read
  [the developer handbook](developer/index.html),
  [adding models](usage/newmodels/index.html), and the
  [contributor guide](contributing/index.html).

## How to read the examples

Complete programs can be copied into a project's `app/main.f90`. A block
explicitly described as a snippet or an insertable block depends on the setup
shown before it; do not compile it as a standalone program. Routine names
linked into the generated API lead to argument-level documentation.

Most Fortran examples use `pr`, the library's real-number precision, and the
default unit system: pressure in bar, volume in litres, temperature in kelvin,
and amount in moles. Check model parameter units separately. Component arrays
must always have the same ordering as the model's construction data.

The guides introduce the physical idea first, then a practical calculation,
then the equations or advanced options. It is fine to postpone the detailed
derivations until you need them. No model or convergence status can replace
checking units, parameters, material balances, and comparison with measured
data.

This library and its documentation are evolving. Please report problems or
missing explanations through the
[issue tracker](https://github.com/ipqa-research/yaeos/issues), including your
library version and a small reproducing example when possible.
