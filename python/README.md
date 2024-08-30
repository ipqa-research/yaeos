# yaeos Python bindings
THIS IS A WIP SO THE API WILL DRASTICALLY CHANGE WITH TIME

Set of Python bindings to call `yaeos` functions and models.

Editable installation

```
cd python
pip install -r requirements-build.txt
pip install -e . --no-build-isolation
```

If you want to install on your environment instead

```shell
pip install .
```

To check if the installation worked correctly:

```python
from yaeos import PengRobinson76

import numpy as np

model = PengRobinson76(
    np.array([320, 375]),
    np.array([30, 45]),
    np.array([0.0123, 0.045])
)

model.lnphi_vt(np.array([5.0, 4.0]), 2.0, 303.15)
```

```
array([0.47647471, 0.35338115])
```