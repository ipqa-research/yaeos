import os

import pandas as pd

import pathlib

import pytest


PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
DATA_PATH = PATH / "datasets"


@pytest.fixture
def data_path():
    return DATA_PATH.joinpath


@pytest.fixture
def data_co2_c6_pxy(data_path):
    tc = [304.1, 504.0]
    pc = [73.75, 30.12]
    w = [0.4, 0.299]
    return pd.read_csv(data_path("co2_c6.csv")), tc, pc, w
