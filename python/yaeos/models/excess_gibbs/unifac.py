"""UNIFAC Module."""

from typing import List

from yaeos.core import GeModel
from yaeos.lib import yaeos_c
from yaeos.models.groups import groups_from_dicts


class UNIFACVLE(GeModel):
    """UNIFAC VLE model.

    Please refer to the `yaeos` user documentation for an in-depth look at the
    model's information: https://ipqa-research.github.io/yaeos/page/index.html

    Parameters
    ----------
    molecules : list of dict
        List of dicts with the groups and their amounts for each molecule.

    Example
    -------
    .. code-block:: python

        from yaeos import UNIFACVLE

        # Groups for water and ethanol
        water = {16: 1}
        ethanol = {1: 1, 2: 1, 14: 1}

        groups = [water, ethanol]

        model = UNIFAVLE(groups)

        model.ln_gamma([0.5, 0.5], 298.15)
    """

    def __init__(self, molecules: List[dict]) -> None:

        self.molecules = molecules

        number_of_groups, groups_ids, groups_ammounts = groups_from_dicts(
            molecules
        )
        self.id = yaeos_c.unifac_vle(
            ngs=number_of_groups, g_ids=groups_ids, g_v=groups_ammounts
        )

    def size(self) -> int:
        """Get the number of components.

        Returns
        -------
        int
            Number of components
        """
        return len(self.molecules)

    def _model_params_as_str(self) -> str:
        """Return the model parameters as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters. This string should be valid
        Fortran code that assigns the model variables.
        """
        fcode = ""

        for i in range(self.size()):
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
        fcode += "ge_model = setup_unifac(molecules)\n\n"

        fcode = fcode.replace(", ]", "]")

        return fcode

    def _model_params_declaration_as_str(self) -> str:
        """Return the model parameters declaration as a string.

        This method should be implemented by subclasses to return a string
        representation of the model parameters declaration. This string should
        be valid Fortran code that declares the model variables.
        """
        fcode = (
            f"integer, parameter :: nc={self.size()}\n"
            "\n"
            "type(UNIFAC) :: ge_model\n"
            "\n"
            f"type(Groups) :: molecules(nc)\n\n"
        )

        return fcode
