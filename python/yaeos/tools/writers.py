"""Module of helper tools that write data to either files or strings in various 
formats.
"""


import numpy as np


def fmatrix_as_str(array: np.ndarray) -> str:
    """Converts a 2D array into it's equivalent representation in Fortran code.
    """

    array_str = ""
    n = len(array)

    for i in range(n):
        array_str += f"array({i + 1}, :) = ["

        for j in range(n):
            if j < n - 1:
                array_str += f"{array[i][j]}_pr, "
            else:
                array_str += f"{array[i][j]}_pr]\n"
    return array_str
