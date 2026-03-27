import subprocess
from pathlib import Path
from typing import List

import numpy as np

import pytest

import yaeos
from yaeos.core import ArModel

# Path to fortran API app folder
here = Path(__file__).parent.parent.parent.parent.parent / "app"

# Define models and add them to the models list
# Models should have 3 compunds
models_list: List[ArModel] = []
# For each model add an identifier to the model_names list
model_names: List[str] = []

# =============================================================================
# Cubic equations of state
# -----------------------------------------------------------------------------
# Normal cubics
tcs = [591.95, 507.6, 514.0]
pcs = [57.86, 30.25, 61.37]
ws = [0.466521, 0.301261, 0.643558]
zcs = [0.211, 0.266, 0.241]

classes = [
    yaeos.PengRobinson76,
    yaeos.PengRobinson78,
    yaeos.SoaveRedlichKwong,
    yaeos.RKPR,
]

for cla in classes:
    if cla is yaeos.RKPR:
        model = cla(
            critical_temperatures=tcs,
            critical_pressures=pcs,
            acentric_factors=ws,
            critical_z=zcs,
        )
    else:
        model = cla(
            critical_temperatures=tcs,
            critical_pressures=pcs,
            acentric_factors=ws,
        )
    models_list.append(model)
    model_names.append(f"{cla.name} with no mixrule")

# Add RKPR with d1 and k
model = yaeos.RKPR(
    critical_temperatures=tcs,
    critical_pressures=pcs,
    acentric_factors=ws,
    critical_z=zcs,
    k=models_list[-1].k,
    delta_1=models_list[-1].delta_1,
)
models_list.append(model)
model_names.append(f"RKPR with d1 and k with no mixrule")

# With mixrules
mixrules = []

# QMR
kij = np.array([[0.0, 0.01, 0.02], [0.01, 0.0, 0.03], [0.02, 0.03, 0.0]])

qmr = yaeos.QMR(kij=kij, lij=kij / 2)

mixrules.append(qmr)

# QMRTD
qmrtd = yaeos.QMRTD(
    kij_0=kij,
    kij_inf=kij * 0.05,
    t_ref=np.array(
        [[0.0, 250.0, 285.0], [250.0, 0.0, 300.0], [285.0, 300.0, 0.0]]
    ),
    lij=kij / 2,
)

mixrules.append(qmrtd)

# MHV
ge = yaeos.UNIFACVLE([{1: 1, 42: 1}, {1: 2, 2: 4}, {1: 1, 2: 1, 14: 1}])
mhv = yaeos.MHV(ge, -0.55)
mixrules.append(mhv)

# HV
hv = yaeos.HV(ge)
mixrules.append(hv)

# HVNRTL
alpha = np.array([[0.0, 0.2, 0.3], [0.1, 0.0, 0.4], [0.2, 0.3, 0.0]])
gji = kij * 0.01
gjit = kij * 0.001

hv_nrtl = yaeos.HVNRTL(alpha, gji, gjit)
mixrules.append(hv_nrtl)

# BUILD MODELS WITH MIXRULES
for i, mixrule in enumerate(mixrules):
    for j, cla in enumerate(classes):
        if cla is yaeos.RKPR:
            model = cla(
                critical_temperatures=tcs,
                critical_pressures=pcs,
                acentric_factors=ws,
                critical_z=zcs,
                mixrule=mixrule,
            )
        else:
            model = cla(
                critical_temperatures=tcs,
                critical_pressures=pcs,
                acentric_factors=ws,
                mixrule=mixrule,
            )

        models_list.append(model)
        model_names.append(f"{cla.name} with {mixrule.name} mixrule")

# =============================================================================
# Multifluid EoS
# -----------------------------------------------------------------------------
model = yaeos.GERG2008(["n-hexane", "decane", "isobutane"])
models_list.append(model)
model_names.append("GERG2008")

# =============================================================================
# SAFT EoS
# -----------------------------------------------------------------------------
model = yaeos.PCSAFT(
    m=[3.0576, 4.6627, 2.2616],
    sigma=[3.7983, 3.8384, 3.7574],
    epsilon_k=[236.77, 243.87, 216.53],
)
models_list.append(model)
model_names.append("PC-SAFT with no mixrule")

model = yaeos.PCSAFT(
    m=[3.0576, 4.6627, 2.2616],
    sigma=[3.7983, 3.8384, 3.7574],
    epsilon_k=[236.77, 243.87, 216.53],
    kij=np.array([[0.0, 0.01, 0.02], [0.01, 0.0, 0.03], [0.02, 0.03, 0.0]]),
)
models_list.append(model)
model_names.append("PC-SAFT with kij")


# ! ===========================================================================
# ! Evaluate all termoprops fortran code
# ! ---------------------------------------------------------------------------
thermoprops_code = """
eval_thermoprops : block 
    real(pr) :: n(nc), T, P, V
    
    ! Residual properties
    real(pr) :: Ar, ArT, ArP, ArT2, ArV, ArV2, Arn(nc), ArTV, ArTn(nc), ArVn(nc), Arn2(nc, nc)
    real(pr) :: Pm, dPdV, dPdT, dPdn(nc)
    real(pr) :: lnPhi(nc), dlnPhidT(nc), dlnPhidn(nc, nc), dlnPhidP(nc)
    real(pr) :: Hr, HrT, HrV, HrP, Hrn(nc)
    real(pr) :: Sr, SrT, SrV, SrP, Srn(nc)
    real(pr) :: Ur, UrT, UrV, UrP, Urn(nc)
    real(pr) :: Gr, GrT, GrV, GrP, Grn(nc)
    real(pr) :: Cv, Cp
    
    ! Excess properties
    real(pr) :: Ge, GeT, GeP, Gen(nc)
    real(pr) :: Ve
    real(pr) :: lngamma(nc), dlngammadP(nc), dlngammadT(nc), dlngammadn(nc, nc)
    real(pr) :: Ae
    real(pr) :: He
    real(pr) :: Ue
    real(pr) :: Se

    integer :: i

    n = [4.0_pr, 10.15_pr, 12.1_pr]
    T = 298.15_pr
    P = 1.01325_pr
    
    call ar_model%volume(n, T, P, V=V, root_type="stable")


    ! =========================================================================
    ! VT residual properties
    ! -------------------------------------------------------------------------
    ! ! Ar_vt
    call ar_model%residual_helmholtz(&
       n, V, T, Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, &
       Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)
    
    print *, "=== FORTRAN OUTPUT START ==="

    ! Scalars
    write(*,'(f0.12)') Ar
    write(*,'(f0.12)') ArV
    write(*,'(f0.12)') ArT
    write(*,'(f0.12)') ArTV
    write(*,'(f0.12)') ArV2
    write(*,'(f0.12)') ArT2

    ! Vectors
    write(*,'(*(f0.12,:,","))') Arn
    write(*,'(*(f0.12,:,","))') ArVn
    write(*,'(*(f0.12,:,","))') ArTn

    ! Matrix by row
    do i=1,nc
        write(*,'(*(f0.12,:,","))') Arn2(i,:)
    end do



    ! ! Pressure
    call ar_model%pressure(n, V, T, P=P, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)

    ! Scalars
    write(*,'(f0.12)') P
    write(*,'(f0.12)') dPdV
    write(*,'(f0.12)') dPdT

    ! Vectors
    write(*,'(*(f0.12,:,","))') dPdn
    
    

    ! ! ln_phi_vt
    call ar_model%lnphi_vt(n=n, V=V, T=T, lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn)
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') lnPhi
    write(*,'(*(f0.12,:,","))') dlnPhidP
    write(*,'(*(f0.12,:,","))') dlnPhidT
    
    ! Matrix by row
    do i=1,nc
        write(*,'(*(f0.12,:,","))') dlnPhidn(i,:)
    end do
    
    
    
    ! ! Hr_vt
    call ar_model%enthalpy_residual_vt(n, V, T, Hr=Hr, HrT=HrT, HrV=HrV, Hrn=Hrn)
    
    ! Scalars
    write(*,'(f0.12)') Hr
    write(*,'(f0.12)') HrT
    write(*,'(f0.12)') HrV
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Hrn
    
    
    
    ! ! Sr_vt
    call ar_model%entropy_residual_vt(n, V, T, Sr=Sr, SrT=SrT, SrV=SrV, Srn=Srn)
    
    ! Scalars
    write(*,'(f0.12)') Sr
    write(*,'(f0.12)') SrT
    write(*,'(f0.12)') SrV
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Srn
    
    
    
    ! ! Ur_vt
    call ar_model%internal_energy_residual_vt(n, V, T, Ur=Ur, UrT=UrT, UrV=UrV, Urn=Urn)
    
    ! Scalars
    write(*,'(f0.12)') Ur
    write(*,'(f0.12)') UrT
    write(*,'(f0.12)') UrV
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Urn
    
    
    
    ! ! Gr_vt
    call ar_model%gibbs_residual_vt(n, V, T, Gr=Gr, GrT=GrT, GrV=GrV, Grn=Grn)
    
    ! Scalars
    write(*,'(f0.12)') Gr
    write(*,'(f0.12)') GrT
    write(*,'(f0.12)') GrV
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Grn
    
    
    
    ! ! Cv_vt
    call ar_model%Cv_residual_vt(n, V, T, Cv=Cv)
    
    ! Scalar
    write(*,'(f0.12)') Cv
    
    
    
    ! ! Cp_vt
    call ar_model%Cp_residual_vt(n, V, T, Cp=Cp)
    
    ! Scalar
    write(*,'(f0.12)') Cp
    
    ! =========================================================================
    ! PT residual properties
    ! -------------------------------------------------------------------------
    ! ! Ar_pt
    call ar_model%helmholtz_residual_pt(n, P, T, root_type="stable", Ar=Ar, &
       ArP=ArP, ArT=ArT, Arn=Arn)
    
    ! Scalars
    write(*,'(f0.12)') Ar
    write(*,'(f0.12)') ArP
    write(*,'(f0.12)') ArT
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Arn
    
    
    
    ! ! ln_phi_pt
    call ar_model%lnphi_pt(n=n, P=P, T=T, root_type="stable", lnPhi=lnPhi, &
      dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn)
      
    ! Vectors
    write(*,'(*(f0.12,:,","))') lnPhi
    write(*,'(*(f0.12,:,","))') dlnPhidP
    write(*,'(*(f0.12,:,","))') dlnPhidT
    
    ! Matrix by row
    do i=1,nc
        write(*,'(*(f0.12,:,","))') dlnPhidn(i,:)
    end do
    
    
    
    ! ! Hr_pt
    call ar_model%enthalpy_residual_pt(n, P, T, root_type="stable", &
       Hr=Hr, HrT=HrT, HrP=HrP, Hrn=Hrn)
    
    ! Scalars
    write(*,'(f0.12)') Hr
    write(*,'(f0.12)') HrT
    write(*,'(f0.12)') HrP
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Hrn
    
    
    
    ! ! Sr_pt
    call ar_model%entropy_residual_pt(n, P, T, root_type="stable", &
       Sr=Sr, SrT=SrT, SrP=SrP, Srn=Srn)
    
    ! Scalars
    write(*,'(f0.12)') Sr
    write(*,'(f0.12)') SrT
    write(*,'(f0.12)') SrP
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Srn
    
    
    
    ! ! Ur_pt
    call ar_model%internal_energy_residual_pt(n, P, T, root_type="stable", &
       Ur=Ur, UrT=UrT, UrP=UrP, Urn=Urn)
    
    ! Scalars
    write(*,'(f0.12)') Ur
    write(*,'(f0.12)') UrT
    write(*,'(f0.12)') UrP
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Urn
    
    
    
    ! ! Gr_pt
    call ar_model%gibbs_residual_pt(n, P, T, root_type="stable", Gr=Gr, &
       GrT=GrT, GrP=GrP, Grn=Grn)
    
    ! Scalars
    write(*,'(f0.12)') Gr
    write(*,'(f0.12)') GrT
    write(*,'(f0.12)') GrP
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Grn
    
    
    
    ! ! Cv_pt
    call ar_model%Cv_residual_pt(n, P, T, root_type="stable", Cv=Cv)
    
    ! Scalar
    write(*,'(f0.12)') Cv
    
    
    
    
    
    ! ! Cp_pt
    call ar_model%Cp_residual_pt(n, P, T, root_type="stable",Cp=Cp)
    
    ! Scalar
    write(*,'(f0.12)') Cp
    
    ! =========================================================================
    ! Excess properties
    ! -------------------------------------------------------------------------
    ! ! Ge
    call ar_model%gibbs_excess(n, P, T, root_type="stable", Ge=Ge, GeT=GeT, &
       GeP=GeP, Gen=Gen)
    
    ! Scalar
    write(*,'(f0.12)') Ge
    write(*,'(f0.12)') GeT
    write(*,'(f0.12)') GeP
    
    ! Vectors
    write(*,'(*(f0.12,:,","))') Gen
    
    
    
    ! ! Ve
    call ar_model%volume_excess(n, P, T, root_type="stable", Ve=Ve)
    
    ! Scalar
    write(*,'(f0.12)') Ve
    
    
    
    ! ! lngamma
    call ar_model%ln_activity_coefficient(n, P, T, root_type="stable", &
       lngamma=lngamma, dlngammadP=dlngammadP, dlngammadT=dlngammadT, &
       dlngammadn=dlngammadn)
      
    ! Vectors
    write(*,'(*(f0.12,:,","))') lngamma
    write(*,'(*(f0.12,:,","))') dlngammadP
    write(*,'(*(f0.12,:,","))') dlngammadT
    
    ! Matrix by row
    do i=1,nc
       write(*,'(*(f0.12,:,","))') dlngammadn(i,:)
    end do
    
    
    
    ! ! Ae
    call ar_model%helmholtz_excess(n, P, T, root_type="stable", Ae=Ae)
    
    ! Scalar
    write(*,'(f0.12)') Ae
    
    
    
    ! ! He
    call ar_model%enthalpy_excess(n, P, T, root_type="stable", He=He)
    
    ! Scalar
    write(*,'(f0.12)') He
    
    
    
    ! ! Ue
    call ar_model%internal_energy_excess(n, P, T, root_type="stable", Ue=Ue)
    
    ! Scalar
    write(*,'(f0.12)') Ue
    
    
    ! ! Se
    call ar_model%entropy_excess(n, P, T, root_type="stable", Se=Se)
    
    ! Scalar
    write(*,'(f0.12)') Se
end block eval_thermoprops
"""  # noqa


@pytest.mark.parametrize("ar_model, name", list(zip(models_list, model_names)))
def test_ar_same_as_fortran(ar_model: ArModel, name: str):
    # Build code
    fortran_code = ""

    fortran_code += "program main\n"
    fortran_code += "use yaeos\n\n"
    fortran_code += "implicit none\n\n"
    fortran_code += ar_model._model_params_declaration_as_str()
    fortran_code += ar_model._model_params_as_str()
    fortran_code += thermoprops_code
    fortran_code += "end program main"

    # Save .f90 code
    file_name = "test_ar_program"

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

    start_idx = None
    for i, line in enumerate(lines):
        if "=== FORTRAN OUTPUT START ===" in line:
            start_idx = i + 1
            break

    # If compilation error occured, here will be printed
    assert (
        start_idx is not None
    ), f"No FORTRAN OUTPUT found.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"  # noqa

    # Actual output lines
    lines = lines[start_idx:]

    # =========================================================================
    # Get values
    # =========================================================================
    # Ar Vt
    ar_f = float(lines[0])
    arv_f = float(lines[1])
    art_f = float(lines[2])
    artv_f = float(lines[3])
    arv2_f = float(lines[4])
    art2_f = float(lines[5])
    arn_f = np.array([float(x) for x in lines[6].split(",")], dtype=float)
    arvn_f = np.array([float(x) for x in lines[7].split(",")], dtype=float)
    artn_f = np.array([float(x) for x in lines[8].split(",")], dtype=float)

    arn2_f = np.array(
        [[float(x) for x in line.split(",")] for line in lines[9:12]],
        dtype=float,
    )

    # =========================================================================
    # Start testing
    # -------------------------------------------------------------------------
    n = np.array([4.0, 10.15, 12.1])
    t = 298.15
    p = 1.01325
    v = ar_model.volume(n, t, p, root="stable")

    # =========================================================================
    # Ar Vt
    # -------------------------------------------------------------------------
    ar, ar_d = ar_model.helmholtz_residual_vt(
        n,
        v,
        t,
        dt=True,
        dv=True,
        dn=True,
        dtv=True,
        dv2=True,
        dt2=True,
        dvn=True,
        dtn=True,
        dn2=True,
    )

    assert np.isclose(ar, ar_f), f"Ar_vt {name}"
    assert np.isclose(ar_d["dt"], art_f), f"ArT_vt {name}"
    assert np.isclose(ar_d["dv"], arv_f), f"ArV_vt {name}"
    assert np.isclose(ar_d["dtv"], artv_f), f"ArTV_vt {name}"
    assert np.isclose(ar_d["dv2"], arv2_f), f"ArV2_vt {name}"
    assert np.isclose(ar_d["dt2"], art2_f), f"ArT2_vt {name}"
    assert np.allclose(ar_d["dn"], arn_f), f"ArN_vt {name}"
    assert np.allclose(ar_d["dvn"], arvn_f), f"ArVN_vt {name}"
    assert np.allclose(ar_d["dtn"], artn_f), f"ArTN_vt {name}"
    assert np.allclose(ar_d["dn2"], arn2_f), f"ArN2_vt {name}"

    # individuals
    # ge_i = ge_model.excess_gibbs(n, t)
    # assert np.isclose(ge_i, ge_f)

    # ge_i, ge_d = ge_model.excess_gibbs(n, t, dt=True)
    # assert np.isclose(ge_i, ge_f)
    # assert np.isclose(ge_d["dt"], get_f)

    # ge_i, ge_d = ge_model.excess_gibbs(n, t, dt2=True)
    # assert np.isclose(ge_i, ge_f)
    # assert np.isclose(ge_d["dt2"], get2_f)

    # ge_i, ge_d = ge_model.excess_gibbs(n, t, dn=True)
    # assert np.isclose(ge_i, ge_f)
    # assert np.allclose(ge_d["dn"], gen_f)

    # ge_i, ge_d = ge_model.excess_gibbs(n, t, dtn=True)
    # assert np.isclose(ge_i, ge_f)
    # assert np.allclose(ge_d["dtn"], getn_f)

    # ge_i, ge_d = ge_model.excess_gibbs(n, t, dn2=True)
    # assert np.isclose(ge_i, ge_f)
    # assert np.allclose(ge_d["dn2"], gen2_f)
