"""UNIFAC Dortmund Module."""

from typing import List

from yaeos.core import GeModel
from yaeos.lib import yaeos_c
from yaeos.models.groups import groups_from_dicts


class UNIFACDortmund(GeModel):
    """UNIFAC Dortmund model.

    Please refer to the `yaeos` user documentation for an in-depth look at the
    model's information: https://ipqa-research.github.io/yaeos/page/index.html

    Parameters
    ----------
    molecules : list of dict
        List of dicts with the groups and their amounts for each molecule.

    Example
    -------
    .. code-block:: python

        from yaeos import UNIFACDortmund

        # Groups for water and ethanol
        water = {16: 1}
        ethanol = {1: 1, 2: 1, 14: 1}

        groups = [water, ethanol]

        model = UNIFACDortmund(groups)

        model.ln_gamma([0.5, 0.5], 298.15)
    """

    def __init__(self, molecules: List[dict]) -> None:

        (number_of_groups, groups_ids, groups_ammounts) = groups_from_dicts(
            molecules
        )
        self.id = yaeos_c.unifac_dortmund(
            ngs=number_of_groups, g_ids=groups_ids, g_v=groups_ammounts
        )

        self.nc = len(molecules)
        self.g_ids = groups_ids
        self.g_ammounts = groups_ammounts

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        for i in range(self.nc):
            id_c = f"molecules({i + 1})%groups_ids = ["
            am_c = f"molecules({i + 1})%number_of_groups = ["

            for j in range(len(self.g_ammounts[i])):
                if self.g_ammounts[i][j] == 0 or self.g_ids[i][j] == 0:
                    continue

                if j < len(self.g_ammounts) - 1:
                    id_c += f"{self.g_ids[i][j]}, "
                    am_c += f"{self.g_ammounts[i][j]}, "
                else:
                    id_c += f"{self.g_ids[i][j]}"
                    am_c += f"{self.g_ammounts[i][j]}"

            id_c += "]"
            am_c += "]"

            fcode += id_c + "\n" + am_c + "\n\n"

        fcode += "\n"
        fcode += "ge_model = setup_dortmund(molecules)\n\n"

        fcode = fcode.replace(", ]", "]")

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        fcode = (
            f"integer, parameter :: nc={self.nc}\n"
            "\n"
            "type(UNIFAC) :: ge_model\n"
            "\n"
            f"type(Groups) :: molecules(nc)\n\n"
        )

        return fcode
