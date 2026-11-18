---
title: Becoming a yaeos developer
ordered_subpage: architecture.md
                 workflow.md
                 fortran.md
                 numerics.md
                 c_interface.md
                 python_interface.md
                 testing.md
                 walkthrough.md
---

[TOC]

# How to become a yaeos developer

This handbook is for someone who can use a model but now wants to understand
the implementation, fix a calculation, add a feature, or expose it to another
language. It covers the Fortran library **and** the interface layers around it.
You do not need to read every source file before making a useful contribution.

If you have never run `yaeos`, first complete the
[user tutorial](../usage/index.html). For terminology such as residual energy,
phase fraction, or volume root, use
[thermodynamics essentials](../usage/fundamentals.html). This handbook explains
where those concepts live in the code and how to change them safely.

## Your first week: a suggested route

| Step | Guide | Practical outcome |
|------|-------|-------------------|
| 1 | [Architecture and repository map](architecture.html) | Locate the implementation behind a public call |
| 2 | [Development workflow](workflow.html) | Build and test a checkout rather than an installed release |
| 3 | [Working in the Fortran core](fortran.html) | Understand types, bindings, derivative contracts, and exports |
| 4 | [Working on numerical algorithms](numerics.html) | Separate model errors from residual, Jacobian, and solver errors |
| 5 | [The C-facing interface layer](c_interface.html) | Understand integer model handles and the actual interoperability boundary |
| 6 | [The Python interface](python_interface.html) | Follow a Python call through the compiled wrapper into Fortran |
| 7 | [Tests and numerical validation](testing.html) | Add a regression and check more than numerical convergence |
| 8 | [An end-to-end change](walkthrough.html) | Plan a feature across implementation, interfaces, tests, and docs |

## Pick the smallest useful change

For a first contribution, choose a change with a clear expected result: a
missing explanation, an incorrect argument description, an additional
regression test, or a small wrapper around an existing calculation. An entirely
new thermodynamic model is possible, but it adds both mathematical and
software-design questions at once.

Write down three things before editing:

1. **The contract:** what inputs, units, shapes, and output behavior should be
   preserved or added?
2. **The evidence:** what test or calculation will distinguish the old behavior
   from the corrected behavior?
3. **The affected layers:** is this a Fortran-only change, or must callers,
   wrappers, Python classes, exports, and documentation also change?

These questions keep a local fix from becoming an accidental API break.

## A feature can cross several layers

```text
Python user code
    -> Python model method / result conversion
    -> compiled f2py wrapper
    -> yaeos_c module: handle lookup and argument adaptation
    -> Fortran model or equilibrium entry point
    -> property/derivative kernel and numerical routines
```

A native Fortran program starts at the Fortran entry point and bypasses the
Python layers. The directory name `c_interface` and use of `iso_c_binding`
do not, on their own, establish that every routine is an independently callable
standard C ABI. The [interface guide](c_interface.html) explains what the
current build actually exposes and what a new C client would require.

## What this handbook treats as authoritative

The explanation follows the current repository's code and build configuration.
The architecture discussion was informed by the author's thesis material on
modularity, derivative adapters, and language interfaces. The thesis was used
as read-only background, not copied into the repository or treated as an API
specification. Names, signatures, and interface mechanisms can evolve between
a research description and an implementation.

Source comments and generated API pages help you navigate, but verify an
important claim against the actual procedure body and its tests. In
particular, do not assume that every solver has the same failure convention,
that a declared optional argument is handled safely on every path, or that a
model with mutable global support state is thread-safe.

## Before asking for review

Keep implementation, regression tests, and documentation together. Explain
which commands you actually ran and which checks were unavailable. Do not
report a Python build as tested merely because the native Fortran build
succeeded. The [testing guide](testing.html) and the
[contribution guide](../contributing/index.html) provide the next steps.
