---
title: Phase Equilibrium
---

Phase Equilibria calculations are fundamental for the majority of EoS based
modelling either for processes or when studying phase-behaviour.

In `yaeos` most of phase-equilibria procedures return the `EquilibriaState`
type [[EquilibriaState]], which holds all the relevant information of an
equilibria point.

The implemented methods, and their usage are:

[TOC]

# Flash calculations
Flash calcuations are one of the most used phase-equilibria calculations during.

In `yaeos` it is possible to make Flash calculations either specifying:

- \(zPT\)
- \(zVT\)

```fortran
type(EquilibriaState) :: result

! zPT flash
result = flash(model, z, p_spec=P, T=T)

! zVT flash
result = flash(model, z, v_spec=P, T=T)

! It is possible to provide initialization compositions in terms of the
! K-factors. Where k0=y/x
result = flash(model, z, v_spec=P, T=T, k0=k0)
```

# Saturation points
Single saturation point calculations are included with the procedures
[[saturation_pressure]] and [[saturation_temperature]]. Both procedures solve
the equation

\[
f(T \lor P) = \sum ln K_i - 1 = 0
\]

With a newton procedure with respect to the desired variable (either \(P\) or
\(T\).

```fortran
type(EquilibriaState) :: sat_point

sat = saturation_pressure(model, z, T=T, kind="bubble")
sat = saturation_pressure(model, z, T=T, kind="dew")

sat = saturation_temperature(model, z, P=P, kind="bubble")
sat = saturation_temperature(model, z, P=P, kind="dew")
```

## Phase envelopes
Phase envelopes are the conection of all the saturation points of a system.
When the interest is in calculating a whole phase diagram instead of a single
point, or the point is hard to converge. It is better to use a robust
mathematical algorithm that eases the calcuation providing an easy-to-converge
point and using its information to initialize a next one and continue along the
whole phase-boundary. This can be done with the procedure [[pt_envelope_2ph]]

```fortran
type(PTEnvel2) :: env

sat = saturation_pressure(model, z, T=150._pr, kind="bubble")
env = pt_envelope_2ph(model, z, sat)
```