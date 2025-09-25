import subprocess
from pathlib import Path
from typing import List

import numpy as np

import pytest

import yaeos
from yaeos.core import GeModel


# Path to fortran API app folder
here = Path(__file__).parent.parent.parent.parent.parent / "app"

# Define models and add them to the models list
# Models should have 3 compunds
models_list: List[GeModel] = []

# =============================================================================
# NRTL
# -----------------------------------------------------------------------------
bij = bij = np.array(
    [
        [0.0, 198.39, -166.14],
        [-289.66, 0.0, -652.55],
        [1620.9, 778.64, 0.0],
    ]
)

aij = bij * 0.1

alphaij = np.array(
    [
        [0.0, 0.3702, 0.40115],
        [0.3702, 0.0, 0.33541],
        [0.40115, 0.33541, 0.0],
    ]
)

nrtl_model = yaeos.NRTL(aij, bij, alphaij)
models_list.append(nrtl_model)

# ! ===========================================================================
# ! Evaluate all termoprops fortran code
# ! ---------------------------------------------------------------------------
thermoprops_code = """
eval_thermoprops : block 
    real(pr) :: n(nc), T
    real(pr) :: Ge, Gen(nc), Gen2(nc, nc), GeT, GeT2, GeTn(nc)
    real(pr) :: He, HeT, Hen(nc)
    real(pr) :: Se, SeT, Sen(nc)
    real(pr) :: lngamma(nc), dlngammadT(nc), dlngammadn(nc, nc)

    integer :: i

    n = [4.0_pr, 10.15_pr, 12.1_pr]
    T = 298.15_pr

    ! ! GE
    call ge_model%excess_gibbs(n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2)

    ! Scalars
    write(*,'(f0.12)') Ge
    write(*,'(f0.12)') GeT
    write(*,'(f0.12)') GeT2

    ! Vectors
    write(*,'(*(f0.12,:,","))') Gen
    write(*,'(*(f0.12,:,","))') GeTn

    ! Matrix by row
    do i=1,nc
        write(*,'(*(f0.12,:,","))') Gen2(i,:)
    end do

    ! ! He
    call ge_model%excess_enthalpy(n, T, He=He, HeT=HeT, Hen=Hen)

    ! Scalars
    write(*,'(f0.12)') He
    write(*,'(f0.12)') HeT
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Hen

    ! ! Se
    call ge_model%excess_entropy(n, T, Se=Se, SeT=SeT, Sen=Sen)

    ! Scalars
    write(*,'(f0.12)') Se
    write(*,'(f0.12)') SeT
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Sen

    ! ! ln_gammas
    call ge_model%ln_activity_coefficient(n, T, lngamma=lngamma, dlngammadT=dlngammadT, dlngammadn=dlngammadn)

    ! Vectors
    write(*,'(*(f0.12,:,","))') lngamma
    write(*,'(*(f0.12,:,","))') dlngammadT

    ! Matrix by row
    do i=1,nc
        write(*,'(*(f0.12,:,","))') dlngammadn(i,:)
    end do
end block eval_thermoprops
"""  # noqa


@pytest.mark.parametrize("ge_model", models_list)
def test_ge_same_as_fortran(ge_model: GeModel):
    # Build code
    fortran_code = ""

    fortran_code += "program main\n"
    fortran_code += "use yaeos\n\n"
    fortran_code += "implicit none\n\n"
    fortran_code += ge_model._model_params_declaration_as_str()
    fortran_code += ge_model._model_params_as_str()
    fortran_code += thermoprops_code
    fortran_code += "end program main"

    # Save .f90 code
    file_name = "test_ge_program"

    with open(here / (file_name + ".f90"), "w") as f:
        f.write(fortran_code)

    # run with fpm and capture prints
    result = subprocess.run(
        ["fpm", "run", file_name, "--profile", "release"],
        cwd=here,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,  # strings instead bytes
    )

    lines = result.stdout.splitlines()[1:-1]

    # Get values
    # Ge
    ge_f = float(lines[0])
    get_f = float(lines[1])
    get2_f = float(lines[2])
    gen_f = np.array([float(x) for x in lines[3].split(",")], dtype=float)
    getn_f = np.array([float(x) for x in lines[4].split(",")], dtype=float)

    gen2_f = np.zeros((3, 3))
    gen2_f[0, :] = np.array(
        [float(x) for x in lines[5].split(",")], dtype=float
    )
    gen2_f[1, :] = np.array(
        [float(x) for x in lines[6].split(",")], dtype=float
    )
    gen2_f[2, :] = np.array(
        [float(x) for x in lines[7].split(",")], dtype=float
    )

    # He
    he_f = float(lines[8])
    het_f = float(lines[9])
    hen_f = np.array([float(x) for x in lines[10].split(",")], dtype=float)

    # Se
    se_f = float(lines[11])
    set_f = float(lines[12])
    sen_f = np.array([float(x) for x in lines[13].split(",")], dtype=float)

    # lngamma
    gam_f = np.array([float(x) for x in lines[14].split(",")], dtype=float)
    gamt_f = np.array([float(x) for x in lines[15].split(",")], dtype=float)

    gamn_f = np.zeros((3, 3))
    gamn_f[0, :] = np.array(
        [float(x) for x in lines[16].split(",")], dtype=float
    )
    gamn_f[1, :] = np.array(
        [float(x) for x in lines[17].split(",")], dtype=float
    )
    gamn_f[2, :] = np.array(
        [float(x) for x in lines[18].split(",")], dtype=float
    )

    # Start testing
    n = np.array([4.0, 10.15, 12.1])
    t = 298.15

    # Ge
    ge, ge_d = ge_model.excess_gibbs(
        n, t, dt=True, dt2=True, dn=True, dtn=True, dn2=True
    )

    assert np.isclose(ge, ge_f)
    assert np.isclose(ge_d["dt"], get_f)
    assert np.isclose(ge_d["dt2"], get2_f)
    assert np.allclose(ge_d["dn"], gen_f)
    assert np.allclose(ge_d["dtn"], getn_f)
    assert np.allclose(ge_d["dn2"], gen2_f)

    # Ge individuals
    ge_i = ge_model.excess_gibbs(n, t)
    assert np.isclose(ge_i, ge_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dt=True)
    assert np.isclose(ge_i, ge_f)
    assert np.isclose(ge_d["dt"], get_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dt2=True)
    assert np.isclose(ge_i, ge_f)
    assert np.isclose(ge_d["dt2"], get2_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dn=True)
    assert np.isclose(ge_i, ge_f)
    # assert np.allclose(ge_d["dn"], getn_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dtn=True)
    assert np.isclose(ge_i, ge_f)
    assert np.allclose(ge_d["dtn"], getn_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dn2=True)
    assert np.isclose(ge_i, ge_f)
    assert np.allclose(ge_d["dn2"], gen2_f)
