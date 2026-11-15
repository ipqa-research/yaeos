Writers module
==============

The writer functions return Fortran assignment text for one-, two-, or
three-dimensional parameter arrays. The one- and two-dimensional writers
append the ``_pr`` precision suffix that the Fortran side expects; the
three-dimensional writer writes plain values.

``f1darray_as_str(array, array_name)``
    A Fortran vector assignment.
``fmatrix_as_str(array, array_name)``
    A Fortran matrix assignment.
``f3dmatrix_as_str(array, array_name)``
    A Fortran three-dimensional array assignment.

``array_name`` is the target Fortran variable name, not a file name, and
``array`` is the data to write. The output is source text, not a complete
program: declare the variable and import the correct precision in the
receiving code.

These helpers format values; they do not convert them. Pressures stay in bar,
energies in bar L, and interaction parameters keep whatever convention and
units they already had. Preserve component order on both sides so the
Fortran model means the same mixture as the Python one.

.. automodule:: yaeos.tools.writers
    :members:
    :undoc-members:
    :show-inheritance:
