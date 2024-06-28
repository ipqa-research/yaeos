Para instalar como editable

```
cd python
pip install -e .
python setup.py build_fortran
```

Para instalar como no editable, directamente

```
pip install .
```

Desde el interprete para checkear que todo este bien:

```python
from yaeos import PengRobinson76, QMR

import numpy as np

mr = QMR(np.zeros((2,2)), np.zeros((2,2)))

model = PengRobinson76(np.array([320, 375]), np.array([30, 45]), np.array([0.0123, 0.045]), mr)

model.fugacity(np.array([5.0, 4.0]), 2.0, 303.15)
```

```
{'ln_phi': array([2.61640775, 2.49331419]), 'dt': None, 'dp': None, 'dn': None}
```