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

# =============================================================================
# UNIFACLV
# -----------------------------------------------------------------------------
groups = [{1: 2}, {1: 1, 2: 1, 14: 1}, {28: 1}]

unifac_model = yaeos.UNIFACVLE(groups)

models_list.append(unifac_model)

# =============================================================================
# UNIFACDortmund
# -----------------------------------------------------------------------------
groups = [{1: 2}, {1: 1, 2: 1, 14: 1}, {28: 1}]

unifac_dortmund_model = yaeos.UNIFACDortmund(groups)

models_list.append(unifac_dortmund_model)

# =============================================================================
# UNIFACPSRK
# -----------------------------------------------------------------------------
groups = [{1: 2}, {1: 1, 2: 1, 14: 1}, {28: 1}]

unifac_psrk_model = yaeos.UNIFACPSRK(groups)

models_list.append(unifac_psrk_model)

# =============================================================================
# UNIQUAC
# -----------------------------------------------------------------------------
a = np.array(
    [
        [0.0, -75.46, -60.15],
        [120.20, 0.0, 44.22],
        [120.20, 33.21, 0.0],
    ]
)

b = np.array(
    [
        [0.0, -0.10062, 0.2566],
        [0.44835, 0.0, -0.01325],
        [0.44835, 0.124, 0.0],
    ]
)

c = np.array(
    [
        [0.0, -0.0008052, 0.00021],
        [0.0004704, 0.0, -0.00033],
        [0.0004704, -0.000247, 0.0],
    ]
)

d = np.array(
    [
        [0.0, -0.001, 0.0002],
        [-0.001, 0.0, 0.0002],
        [-0.001, 0.0002, 0.0],
    ]
)

e = np.array(
    [
        [0.0, -0.00001, 0.00001],
        [-0.00001, 0.0, 0.00001],
        [-0.00001, 0.00001, 0.0],
    ]
)

rs = np.array([0.92, 2.1055, 1.5])
qs = np.array([1.4, 1.972, 1.4])

uniquac_model = yaeos.UNIQUAC(qs, rs, a, b, c, d, e)

models_list.append(uniquac_model)


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
    real(pr) :: Cpe

    integer :: i

    n = [4.0_pr, 10.15_pr, 12.1_pr]
    T = 298.15_pr

    ! ! GE
    call ge_model%excess_gibbs(n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2)
    
    print *, "=== FORTRAN OUTPUT START ==="

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
    
    ! ! Cpe
    call ge_model%excess_Cp(n, T, Cpe=Cpe)

    ! Scalars
    write(*,'(f0.12)') Cpe
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

    # Search output start
    lines = result.stdout.splitlines()

    start_idx = 0
    for i, line in enumerate(lines):
        if "=== FORTRAN OUTPUT START ===" in line:
            start_idx = i + 1
            break

    # Actual output lines
    lines = lines[start_idx:]

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

    # Cpe
    cpe_f = float(lines[19])

    # Start testing
    n = np.array([4.0, 10.15, 12.1])
    t = 298.15

    # =========================================================================
    # Ge
    # -------------------------------------------------------------------------
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
    assert np.allclose(ge_d["dn"], gen_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dtn=True)
    assert np.isclose(ge_i, ge_f)
    assert np.allclose(ge_d["dtn"], getn_f)

    ge_i, ge_d = ge_model.excess_gibbs(n, t, dn2=True)
    assert np.isclose(ge_i, ge_f)
    assert np.allclose(ge_d["dn2"], gen2_f)

    # =========================================================================
    # He
    # -------------------------------------------------------------------------
    he, he_d = ge_model.excess_enthalpy(n, t, dt=True, dn=True)

    assert np.isclose(he, he_f)
    assert np.isclose(he_d["dt"], het_f)
    assert np.allclose(he_d["dn"], hen_f)

    # He individuals
    he_i = ge_model.excess_enthalpy(n, t)
    assert np.isclose(he_i, he_f)

    he_i, he_d = ge_model.excess_enthalpy(n, t, dt=True)
    assert np.isclose(he_i, he_f)
    assert np.isclose(he_d["dt"], het_f)

    he_i, he_d = ge_model.excess_enthalpy(n, t, dn=True)
    assert np.isclose(he_i, he_f)
    assert np.allclose(he_d["dn"], hen_f)

    # =========================================================================
    # Se
    # -------------------------------------------------------------------------
    se, se_d = ge_model.excess_entropy(n, t, dt=True, dn=True)

    assert np.isclose(se, se_f)
    assert np.isclose(se_d["dt"], set_f)
    assert np.allclose(se_d["dn"], sen_f)

    # se individuals
    se_i = ge_model.excess_entropy(n, t)
    assert np.isclose(se_i, se_f)

    se_i, se_d = ge_model.excess_entropy(n, t, dt=True)
    assert np.isclose(se_i, se_f)
    assert np.isclose(se_d["dt"], set_f)

    se_i, se_d = ge_model.excess_entropy(n, t, dn=True)
    assert np.isclose(se_i, se_f)
    assert np.allclose(se_d["dn"], sen_f)

    # =========================================================================
    # gamma
    # -------------------------------------------------------------------------
    gamma, gamma_d = ge_model.ln_gamma(n, t, dt=True, dn=True)

    assert np.allclose(gamma, gam_f)
    assert np.allclose(gamma_d["dt"], gamt_f)
    assert np.allclose(gamma_d["dn"], gamn_f)

    # gamma individuals
    gamma_i = ge_model.ln_gamma(n, t)
    assert np.allclose(gamma_i, gam_f)

    gamma_i, gamma_d = ge_model.ln_gamma(n, t, dt=True)
    assert np.allclose(gamma_i, gam_f)
    assert np.allclose(gamma_d["dt"], gamt_f)

    gamma_i, gamma_d = ge_model.ln_gamma(n, t, dn=True)
    assert np.allclose(gamma_i, gam_f)
    assert np.allclose(gamma_d["dn"], gamn_f)

    # =========================================================================
    # Cpe
    # -------------------------------------------------------------------------
    cpe = ge_model.excess_cp(n, t)
    assert np.isclose(cpe, cpe_f)
