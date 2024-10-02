# yaeos Python bindings

The `yaeos` Python API is on and early stage of development. So the API will 
change with time.

## Supported operative systems

- Linux

## Supported Python versions

- Python >= 3.10, < 3.13

## Installation
To install the last version of `yaeos` Python API, you can use `pip`:

```
pip install yaeos
```

### For developers, the editable installation is done by

```
cd python
pip install -r requirements-build.txt
pip install -e . --no-build-isolation
```

### Building from source

To build `yaeos` from source you can clone the repository and install it with 
pip. The system dependencies are:

- LAPACK
- Fortran compiler

To build and install do:

```shell
git clone https://github.com/ipqa-research/yaeos.git
cd yaeos/python
pip install .
```

### Setting a Google Colab to use yaeos

To install `yaeos` on Google Colab you can do the following command to install 
the last version uploaded to PyPI:

```shell
%pip install yaeos
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

You will get:

```
array([0.47647471, 0.35338115])
```