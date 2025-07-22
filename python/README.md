# yaeos Python API

![logo](https://github.com/ipqa-research/yaeos/blob/main/media/logo.png?raw=true)

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

For developers, editable installation can be done with:

```shell
git clone https://github.com/ipqa-research/yaeos.git

cd yaeos/python
pip install -r requirements-build.txt
pip install -e . --no-build-isolation
```

### Setting a Google Colab to use yaeos

To install `yaeos` on Google Colab you can do the following command to install 
the last version uploaded to PyPI:

```shell
%pip install yaeos
```
