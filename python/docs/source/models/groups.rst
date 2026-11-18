Groups module
=============

Group-contribution models (the UNIFAC family and PSRK) represent a molecule
as counts of functional subgroups. Python constructors accept one dictionary
per component: keys are integer subgroup IDs and values are counts within
that molecule, **not mixture mole fractions**. The list order establishes the
component order for later property calls.

``groups_from_dicts`` converts these dictionaries to the low-level array
layout; ordinary users can pass dictionaries directly to the model instead.
Use the subgroup definitions and interaction tables for the selected UNIFAC
variant. A chemically plausible decomposition with the wrong IDs can still
describe the wrong molecule.

See the :doc:`UNIFAC tutorial </tutorial/gibbs_excess/models/unifac_models>`
before preparing your own dictionaries.

.. automodule:: yaeos.models.groups
    :members:
    :undoc-members:
    :show-inheritance:
