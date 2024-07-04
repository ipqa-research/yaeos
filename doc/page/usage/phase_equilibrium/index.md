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
Flash calcuations

# Saturation points
Single saturation point calculations are included with the procedures
[[saturation_pressure]] and [[saturation_temperature]]. Both procedures solve
the equation

\[
f(T \lor P) = \sum ln K_i - 1 = 0
\]

With a newton procedure with respect to the desired variable (either \(P\) or
\(T\).

## Phase envelopes
Phase envelopes are the conection of all the saturation points of a system.
When the interest is in calculating a whole phase diagram instead of a single
point, or the point is hard to converge. It is better to use a robust
mathematical algorithm that eases the calcuation providing an easy-to-converge
point and using its information to initialize a next one and continue along the
whole phase-boundary. This can be done with the procedure [[pt_envelope_2ph]]