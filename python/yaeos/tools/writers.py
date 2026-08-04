"""Module tool: write data to either files or strings in various formats."""

import numpy as np


def fmatrix_as_str(array: np.ndarray, array_name: str) -> str:
    """Convert 2D array into its equivalent representation in Fortran code.

    Parameters
    ----------
    array: 2D array
        2D array with the values to write as a Fortran matrix with yaeos
        precision variable `pr`
    array_name: str
        Name that the variable should have in the Fortran side
    """
    array_str = ""
    n = len(array)

    for i in range(n):
        array_str += f"{array_name}({i + 1}, :) = ["

        for j in range(n):
            if j < n - 1:
                array_str += f"{array[i][j]}_pr, "
            else:
                array_str += f"{array[i][j]}_pr]\n"
    return array_str
