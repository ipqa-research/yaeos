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
model_names.append("RKPR with d1 and k with no mixrule")

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
    m=[5.76251925, 3.0322045, 2.77409],
    sigma=[2.29910178, 2.00975607, 3.2557],
    epsilon_k=[189.68460534, 311.87818938, 253.406],
)
models_list.append(model)
model_names.append("PC-SAFT with no mixrule")

model = yaeos.PCSAFT(
    m=[5.76251925, 3.0322045, 2.77409],
    sigma=[2.29910178, 2.00975607, 3.2557],
    epsilon_k=[189.68460534, 311.87818938, 253.406],
    kij=np.array(
        [
            [0.0, -0.00500098, 0.002],
            [-0.00500098, 0.0, 0.003],
            [0.002, 0.003, 0.0],
        ]
    ),
)
models_list.append(model)
model_names.append("PC-SAFT with kij")


# ! ===========================================================================
# ! Evaluate all termoprops fortran code
# ! ---------------------------------------------------------------------------
thermoprops_code = """
eval_thermoprops : block 
    real(pr) :: n(nc), T, P, V, Pcal
    
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
    
    call ar_model%volume(n, P, T, V=V, root_type="stable")


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
    call ar_model%pressure(n, V, T, P=Pcal, dPdV=dPdV, dPdT=dPdT, dPdn=dPdn)

    ! Scalars
    write(*,'(f0.12)') Pcal
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


def parse_vec(line):
    return np.array([float(x) for x in line.split(",")], dtype=float)


def parse_mat(lines):
    return np.array(
        [[float(x) for x in li.split(",")] for li in lines], dtype=float
    )


@pytest.mark.parametrize("ar_model, name", list(zip(models_list, model_names)))
def test_ar_same_as_fortran(ar_model: ArModel, name: str):
    # if name != "PC-SAFT with no mixrule":
    #     pytest.skip()

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
    nc = 3
    idx = 0

    # =========================================================================
    # Parse Fortran output
    # =========================================================================

    # --- Ar VT ---
    ar_f = float(lines[idx])
    idx += 1
    arv_f = float(lines[idx])
    idx += 1
    art_f = float(lines[idx])
    idx += 1
    artv_f = float(lines[idx])
    idx += 1
    arv2_f = float(lines[idx])
    idx += 1
    art2_f = float(lines[idx])
    idx += 1
    arn_f = parse_vec(lines[idx])
    idx += 1
    arvn_f = parse_vec(lines[idx])
    idx += 1
    artn_f = parse_vec(lines[idx])
    idx += 1
    arn2_f = parse_mat(lines[idx : idx + nc])
    idx += nc

    # --- Pressure VT ---
    p_f = float(lines[idx])
    idx += 1
    dpdv_f = float(lines[idx])
    idx += 1
    dpdt_f = float(lines[idx])
    idx += 1
    dpdn_f = parse_vec(lines[idx])
    idx += 1

    # --- lnPhi VT ---
    lnphi_vt_f = parse_vec(lines[idx])
    idx += 1
    dlnphidp_vt_f = parse_vec(lines[idx])
    idx += 1
    dlnphidt_vt_f = parse_vec(lines[idx])
    idx += 1
    dlnphidn_vt_f = parse_mat(lines[idx : idx + nc])
    idx += nc

    # --- Hr VT ---
    hr_f = float(lines[idx])
    idx += 1
    hrt_f = float(lines[idx])
    idx += 1
    hrv_f = float(lines[idx])
    idx += 1
    hrn_f = parse_vec(lines[idx])
    idx += 1

    # --- Sr VT ---
    sr_f = float(lines[idx])
    idx += 1
    srt_f = float(lines[idx])
    idx += 1
    srv_f = float(lines[idx])
    idx += 1
    srn_f = parse_vec(lines[idx])
    idx += 1

    # --- Ur VT ---
    ur_f = float(lines[idx])
    idx += 1
    urt_f = float(lines[idx])
    idx += 1
    urv_f = float(lines[idx])
    idx += 1
    urn_f = parse_vec(lines[idx])
    idx += 1

    # --- Gr VT ---
    gr_f = float(lines[idx])
    idx += 1
    grt_f = float(lines[idx])
    idx += 1
    grv_f = float(lines[idx])
    idx += 1
    grn_f = parse_vec(lines[idx])
    idx += 1

    # --- Cv VT ---
    cv_vt_f = float(lines[idx])
    idx += 1

    # --- Cp VT ---
    cp_vt_f = float(lines[idx])
    idx += 1

    # --- Ar PT ---
    ar_pt_f = float(lines[idx])
    idx += 1
    arp_f = float(lines[idx])
    idx += 1
    art_pt_f = float(lines[idx])
    idx += 1
    arn_pt_f = parse_vec(lines[idx])
    idx += 1

    # --- lnPhi PT ---
    lnphi_pt_f = parse_vec(lines[idx])
    idx += 1
    dlnphidp_pt_f = parse_vec(lines[idx])
    idx += 1
    dlnphidt_pt_f = parse_vec(lines[idx])
    idx += 1
    dlnphidn_pt_f = parse_mat(lines[idx : idx + nc])
    idx += nc

    # --- Hr PT ---
    hr_pt_f = float(lines[idx])
    idx += 1
    hrt_pt_f = float(lines[idx])
    idx += 1
    hrp_f = float(lines[idx])
    idx += 1
    hrn_pt_f = parse_vec(lines[idx])
    idx += 1

    # --- Sr PT ---
    sr_pt_f = float(lines[idx])
    idx += 1
    srt_pt_f = float(lines[idx])
    idx += 1
    srp_f = float(lines[idx])
    idx += 1
    srn_pt_f = parse_vec(lines[idx])
    idx += 1

    # --- Ur PT ---
    ur_pt_f = float(lines[idx])
    idx += 1
    urt_pt_f = float(lines[idx])
    idx += 1
    urp_f = float(lines[idx])
    idx += 1
    urn_pt_f = parse_vec(lines[idx])
    idx += 1

    # --- Gr PT ---
    gr_pt_f = float(lines[idx])
    idx += 1
    grt_pt_f = float(lines[idx])
    idx += 1
    grp_f = float(lines[idx])
    idx += 1
    grn_pt_f = parse_vec(lines[idx])
    idx += 1

    # --- Cv PT ---
    cv_pt_f = float(lines[idx])
    idx += 1

    # --- Cp PT ---
    cp_pt_f = float(lines[idx])
    idx += 1

    # --- Ge ---
    ge_f = float(lines[idx])
    idx += 1
    get_f = float(lines[idx])
    idx += 1
    gep_f = float(lines[idx])
    idx += 1
    gen_f = parse_vec(lines[idx])
    idx += 1

    # --- Ve ---
    ve_f = float(lines[idx])
    idx += 1

    # --- lngamma ---
    lngamma_f = parse_vec(lines[idx])
    idx += 1
    dlngammadp_f = parse_vec(lines[idx])
    idx += 1
    dlngammadt_f = parse_vec(lines[idx])
    idx += 1
    dlngammadn_f = parse_mat(lines[idx : idx + nc])
    idx += nc

    # --- Ae ---
    ae_f = float(lines[idx])
    idx += 1

    # --- He ---
    he_f = float(lines[idx])
    idx += 1

    # --- Ue ---
    ue_f = float(lines[idx])
    idx += 1

    # --- Se ---
    se_f = float(lines[idx])
    idx += 1

    # =========================================================================
    # Python evaluations
    # =========================================================================
    n = np.array([4.0, 10.15, 12.1])
    t = 298.15
    p = 1.01325
    v = ar_model.volume(n, p, t, root="stable")

    # =========================================================================
    # Ar VT — vs Fortran + individual derivatives
    # =========================================================================
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

    # Individual derivatives
    ar_i = ar_model.helmholtz_residual_vt(n, v, t)
    assert np.isclose(ar_i, ar_f), f"Ar_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dt=True)
    assert np.isclose(d["dt"], art_f), f"ArT_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dv=True)
    assert np.isclose(d["dv"], arv_f), f"ArV_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dn=True)
    assert np.allclose(d["dn"], arn_f), f"ArN_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dtv=True)
    assert np.isclose(d["dtv"], artv_f), f"ArTV_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dv2=True)
    assert np.isclose(d["dv2"], arv2_f), f"ArV2_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dt2=True)
    assert np.isclose(d["dt2"], art2_f), f"ArT2_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dvn=True)
    assert np.allclose(d["dvn"], arvn_f), f"ArVN_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dtn=True)
    assert np.allclose(d["dtn"], artn_f), f"ArTN_vt individual {name}"

    _, d = ar_model.helmholtz_residual_vt(n, v, t, dn2=True)
    assert np.allclose(d["dn2"], arn2_f), f"ArN2_vt individual {name}"

    # =========================================================================
    # Pressure — vs Fortran + individual derivatives
    # =========================================================================
    p_py, p_d = ar_model.pressure(n, v, t, dv=True, dt=True, dn=True)

    assert np.isclose(p_py, p_f), f"P {name}"
    assert np.isclose(p_d["dv"], dpdv_f), f"dPdV {name}"
    assert np.isclose(p_d["dt"], dpdt_f), f"dPdT {name}"
    assert np.allclose(p_d["dn"], dpdn_f), f"dPdn {name}"

    # Individual
    p_i = ar_model.pressure(n, v, t)
    assert np.isclose(p_i, p_f), f"P individual {name}"

    _, d = ar_model.pressure(n, v, t, dv=True)
    assert np.isclose(d["dv"], dpdv_f), f"dPdV individual {name}"

    _, d = ar_model.pressure(n, v, t, dt=True)
    assert np.isclose(d["dt"], dpdt_f), f"dPdT individual {name}"

    _, d = ar_model.pressure(n, v, t, dn=True)
    assert np.allclose(d["dn"], dpdn_f), f"dPdn individual {name}"

    # =========================================================================
    # lnPhi VT — vs Fortran + individual derivatives
    # =========================================================================
    lnphi_vt, lnphi_vt_d = ar_model.lnphi_vt(
        n, v, t, dt=True, dp=True, dn=True
    )

    assert np.allclose(lnphi_vt, lnphi_vt_f), f"lnPhi_vt {name}"
    assert np.allclose(lnphi_vt_d["dt"], dlnphidt_vt_f), f"dlnPhidT_vt {name}"
    assert np.allclose(lnphi_vt_d["dp"], dlnphidp_vt_f), f"dlnPhidP_vt {name}"
    assert np.allclose(lnphi_vt_d["dn"], dlnphidn_vt_f), f"dlnPhidn_vt {name}"

    # Individual
    lnphi_i = ar_model.lnphi_vt(n, v, t)
    assert np.allclose(lnphi_i, lnphi_vt_f), f"lnPhi_vt individual {name}"

    _, d = ar_model.lnphi_vt(n, v, t, dt=True)
    assert np.allclose(
        d["dt"], dlnphidt_vt_f
    ), f"dlnPhidT_vt individual {name}"

    _, d = ar_model.lnphi_vt(n, v, t, dp=True)
    assert np.allclose(
        d["dp"], dlnphidp_vt_f
    ), f"dlnPhidP_vt individual {name}"

    _, d = ar_model.lnphi_vt(n, v, t, dn=True)
    assert np.allclose(
        d["dn"], dlnphidn_vt_f
    ), f"dlnPhidn_vt individual {name}"

    # =========================================================================
    # Hr VT — vs Fortran + individual derivatives
    # =========================================================================
    hr, hr_d = ar_model.enthalpy_residual_vt(
        n, v, t, dt=True, dv=True, dn=True
    )

    assert np.isclose(hr, hr_f), f"Hr_vt {name}"
    assert np.isclose(hr_d["dt"], hrt_f), f"HrT_vt {name}"
    assert np.isclose(hr_d["dv"], hrv_f), f"HrV_vt {name}"
    assert np.allclose(hr_d["dn"], hrn_f), f"Hrn_vt {name}"

    # Individual
    hr_i = ar_model.enthalpy_residual_vt(n, v, t)
    assert np.isclose(hr_i, hr_f), f"Hr_vt individual {name}"

    _, d = ar_model.enthalpy_residual_vt(n, v, t, dt=True)
    assert np.isclose(d["dt"], hrt_f), f"HrT_vt individual {name}"

    _, d = ar_model.enthalpy_residual_vt(n, v, t, dv=True)
    assert np.isclose(d["dv"], hrv_f), f"HrV_vt individual {name}"

    _, d = ar_model.enthalpy_residual_vt(n, v, t, dn=True)
    assert np.allclose(d["dn"], hrn_f), f"Hrn_vt individual {name}"

    # =========================================================================
    # Sr VT — vs Fortran + individual derivatives
    # =========================================================================
    sr, sr_d = ar_model.entropy_residual_vt(n, v, t, dt=True, dv=True, dn=True)

    assert np.isclose(sr, sr_f), f"Sr_vt {name}"
    assert np.isclose(sr_d["dt"], srt_f), f"SrT_vt {name}"
    assert np.isclose(sr_d["dv"], srv_f), f"SrV_vt {name}"
    assert np.allclose(sr_d["dn"], srn_f), f"Srn_vt {name}"

    # Individual
    sr_i = ar_model.entropy_residual_vt(n, v, t)
    assert np.isclose(sr_i, sr_f), f"Sr_vt individual {name}"

    _, d = ar_model.entropy_residual_vt(n, v, t, dt=True)
    assert np.isclose(d["dt"], srt_f), f"SrT_vt individual {name}"

    _, d = ar_model.entropy_residual_vt(n, v, t, dv=True)
    assert np.isclose(d["dv"], srv_f), f"SrV_vt individual {name}"

    _, d = ar_model.entropy_residual_vt(n, v, t, dn=True)
    assert np.allclose(d["dn"], srn_f), f"Srn_vt individual {name}"

    # =========================================================================
    # Ur VT — vs Fortran + individual derivatives
    # =========================================================================
    ur, ur_d = ar_model.internal_energy_residual_vt(
        n, v, t, dt=True, dv=True, dn=True
    )

    assert np.isclose(ur, ur_f), f"Ur_vt {name}"
    assert np.isclose(ur_d["dt"], urt_f), f"UrT_vt {name}"
    assert np.isclose(ur_d["dv"], urv_f), f"UrV_vt {name}"
    assert np.allclose(ur_d["dn"], urn_f), f"Urn_vt {name}"

    # Individual
    ur_i = ar_model.internal_energy_residual_vt(n, v, t)
    assert np.isclose(ur_i, ur_f), f"Ur_vt individual {name}"

    _, d = ar_model.internal_energy_residual_vt(n, v, t, dt=True)
    assert np.isclose(d["dt"], urt_f), f"UrT_vt individual {name}"

    _, d = ar_model.internal_energy_residual_vt(n, v, t, dv=True)
    assert np.isclose(d["dv"], urv_f), f"UrV_vt individual {name}"

    _, d = ar_model.internal_energy_residual_vt(n, v, t, dn=True)
    assert np.allclose(d["dn"], urn_f), f"Urn_vt individual {name}"

    # =========================================================================
    # Gr VT — vs Fortran + individual derivatives
    # =========================================================================
    gr, gr_d = ar_model.gibbs_residual_vt(n, v, t, dt=True, dv=True, dn=True)

    assert np.isclose(gr, gr_f), f"Gr_vt {name}"
    assert np.isclose(gr_d["dt"], grt_f), f"GrT_vt {name}"
    assert np.isclose(gr_d["dv"], grv_f), f"GrV_vt {name}"
    assert np.allclose(gr_d["dn"], grn_f), f"Grn_vt {name}"

    # Individual
    gr_i = ar_model.gibbs_residual_vt(n, v, t)
    assert np.isclose(gr_i, gr_f), f"Gr_vt individual {name}"

    _, d = ar_model.gibbs_residual_vt(n, v, t, dt=True)
    assert np.isclose(d["dt"], grt_f), f"GrT_vt individual {name}"

    _, d = ar_model.gibbs_residual_vt(n, v, t, dv=True)
    assert np.isclose(d["dv"], grv_f), f"GrV_vt individual {name}"

    _, d = ar_model.gibbs_residual_vt(n, v, t, dn=True)
    assert np.allclose(d["dn"], grn_f), f"Grn_vt individual {name}"

    # =========================================================================
    # Cv / Cp VT — vs Fortran (scalars only, no derivatives)
    # =========================================================================
    cv_vt = ar_model.cv_residual_vt(n, v, t)
    assert np.isclose(cv_vt, cv_vt_f), f"Cv_vt {name}"

    cp_vt = ar_model.cp_residual_vt(n, v, t)
    assert np.isclose(cp_vt, cp_vt_f), f"Cp_vt {name}"

    # =========================================================================
    # Ar PT — vs Fortran + individual derivatives
    # =========================================================================
    ar_pt, ar_pt_d = ar_model.helmholtz_residual_pt(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.isclose(ar_pt, ar_pt_f), f"Ar_pt {name}"
    assert np.isclose(ar_pt_d["dt"], art_pt_f), f"ArT_pt {name}"
    assert np.isclose(ar_pt_d["dp"], arp_f), f"ArP_pt {name}"
    assert np.allclose(ar_pt_d["dn"], arn_pt_f), f"Arn_pt {name}"

    # Individual
    ar_pt_i = ar_model.helmholtz_residual_pt(n, p, t, root="stable")
    assert np.isclose(ar_pt_i, ar_pt_f), f"Ar_pt individual {name}"

    _, d = ar_model.helmholtz_residual_pt(n, p, t, root="stable", dt=True)
    assert np.isclose(d["dt"], art_pt_f), f"ArT_pt individual {name}"

    _, d = ar_model.helmholtz_residual_pt(n, p, t, root="stable", dp=True)
    assert np.isclose(d["dp"], arp_f), f"ArP_pt individual {name}"

    _, d = ar_model.helmholtz_residual_pt(n, p, t, root="stable", dn=True)
    assert np.allclose(d["dn"], arn_pt_f), f"Arn_pt individual {name}"

    # =========================================================================
    # lnPhi PT — vs Fortran + individual derivatives
    # =========================================================================
    lnphi_pt, lnphi_pt_d = ar_model.lnphi_pt(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.allclose(lnphi_pt, lnphi_pt_f), f"lnPhi_pt {name}"
    assert np.allclose(lnphi_pt_d["dt"], dlnphidt_pt_f), f"dlnPhidT_pt {name}"
    assert np.allclose(lnphi_pt_d["dp"], dlnphidp_pt_f), f"dlnPhidP_pt {name}"
    assert np.allclose(lnphi_pt_d["dn"], dlnphidn_pt_f), f"dlnPhidn_pt {name}"

    # Individual
    lnphi_pt_i = ar_model.lnphi_pt(n, p, t, root="stable")
    assert np.allclose(lnphi_pt_i, lnphi_pt_f), f"lnPhi_pt individual {name}"

    _, d = ar_model.lnphi_pt(n, p, t, root="stable", dt=True)
    assert np.allclose(
        d["dt"], dlnphidt_pt_f
    ), f"dlnPhidT_pt individual {name}"

    _, d = ar_model.lnphi_pt(n, p, t, root="stable", dp=True)
    assert np.allclose(
        d["dp"], dlnphidp_pt_f
    ), f"dlnPhidP_pt individual {name}"

    _, d = ar_model.lnphi_pt(n, p, t, root="stable", dn=True)
    assert np.allclose(
        d["dn"], dlnphidn_pt_f
    ), f"dlnPhidn_pt individual {name}"

    # =========================================================================
    # Hr PT — vs Fortran + individual derivatives
    # =========================================================================
    hr_pt, hr_pt_d = ar_model.enthalpy_residual_pt(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.isclose(hr_pt, hr_pt_f), f"Hr_pt {name}"
    assert np.isclose(hr_pt_d["dt"], hrt_pt_f), f"HrT_pt {name}"
    assert np.isclose(hr_pt_d["dp"], hrp_f), f"HrP_pt {name}"
    assert np.allclose(hr_pt_d["dn"], hrn_pt_f), f"Hrn_pt {name}"

    # Individual
    hr_pt_i = ar_model.enthalpy_residual_pt(n, p, t, root="stable")
    assert np.isclose(hr_pt_i, hr_pt_f), f"Hr_pt individual {name}"

    _, d = ar_model.enthalpy_residual_pt(n, p, t, root="stable", dt=True)
    assert np.isclose(d["dt"], hrt_pt_f), f"HrT_pt individual {name}"

    _, d = ar_model.enthalpy_residual_pt(n, p, t, root="stable", dp=True)
    assert np.isclose(d["dp"], hrp_f), f"HrP_pt individual {name}"

    _, d = ar_model.enthalpy_residual_pt(n, p, t, root="stable", dn=True)
    assert np.allclose(d["dn"], hrn_pt_f), f"Hrn_pt individual {name}"

    # =========================================================================
    # Sr PT — vs Fortran + individual derivatives
    # =========================================================================
    sr_pt, sr_pt_d = ar_model.entropy_residual_pt(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.isclose(sr_pt, sr_pt_f), f"Sr_pt {name}"
    assert np.isclose(sr_pt_d["dt"], srt_pt_f), f"SrT_pt {name}"
    assert np.isclose(sr_pt_d["dp"], srp_f), f"SrP_pt {name}"
    assert np.allclose(sr_pt_d["dn"], srn_pt_f), f"Srn_pt {name}"

    # Individual
    sr_pt_i = ar_model.entropy_residual_pt(n, p, t, root="stable")
    assert np.isclose(sr_pt_i, sr_pt_f), f"Sr_pt individual {name}"

    _, d = ar_model.entropy_residual_pt(n, p, t, root="stable", dt=True)
    assert np.isclose(d["dt"], srt_pt_f), f"SrT_pt individual {name}"

    _, d = ar_model.entropy_residual_pt(n, p, t, root="stable", dp=True)
    assert np.isclose(d["dp"], srp_f), f"SrP_pt individual {name}"

    _, d = ar_model.entropy_residual_pt(n, p, t, root="stable", dn=True)
    assert np.allclose(d["dn"], srn_pt_f), f"Srn_pt individual {name}"

    # =========================================================================
    # Ur PT — vs Fortran + individual derivatives
    # =========================================================================
    ur_pt, ur_pt_d = ar_model.internal_energy_residual_pt(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.isclose(ur_pt, ur_pt_f), f"Ur_pt {name}"
    assert np.isclose(ur_pt_d["dt"], urt_pt_f), f"UrT_pt {name}"
    assert np.isclose(ur_pt_d["dp"], urp_f), f"UrP_pt {name}"
    assert np.allclose(ur_pt_d["dn"], urn_pt_f), f"Urn_pt {name}"

    # Individual
    ur_pt_i = ar_model.internal_energy_residual_pt(n, p, t, root="stable")
    assert np.isclose(ur_pt_i, ur_pt_f), f"Ur_pt individual {name}"

    _, d = ar_model.internal_energy_residual_pt(
        n, p, t, root="stable", dt=True
    )
    assert np.isclose(d["dt"], urt_pt_f), f"UrT_pt individual {name}"

    _, d = ar_model.internal_energy_residual_pt(
        n, p, t, root="stable", dp=True
    )
    assert np.isclose(d["dp"], urp_f), f"UrP_pt individual {name}"

    _, d = ar_model.internal_energy_residual_pt(
        n, p, t, root="stable", dn=True
    )
    assert np.allclose(d["dn"], urn_pt_f), f"Urn_pt individual {name}"

    # =========================================================================
    # Gr PT — vs Fortran + individual derivatives
    # =========================================================================
    gr_pt, gr_pt_d = ar_model.gibbs_residual_pt(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.isclose(gr_pt, gr_pt_f), f"Gr_pt {name}"
    assert np.isclose(gr_pt_d["dt"], grt_pt_f), f"GrT_pt {name}"
    assert np.isclose(gr_pt_d["dp"], grp_f), f"GrP_pt {name}"
    assert np.allclose(gr_pt_d["dn"], grn_pt_f), f"Grn_pt {name}"

    # Individual
    gr_pt_i = ar_model.gibbs_residual_pt(n, p, t, root="stable")
    assert np.isclose(gr_pt_i, gr_pt_f), f"Gr_pt individual {name}"

    _, d = ar_model.gibbs_residual_pt(n, p, t, root="stable", dt=True)
    assert np.isclose(d["dt"], grt_pt_f), f"GrT_pt individual {name}"

    _, d = ar_model.gibbs_residual_pt(n, p, t, root="stable", dp=True)
    assert np.isclose(d["dp"], grp_f), f"GrP_pt individual {name}"

    _, d = ar_model.gibbs_residual_pt(n, p, t, root="stable", dn=True)
    assert np.allclose(d["dn"], grn_pt_f), f"Grn_pt individual {name}"

    # =========================================================================
    # Cv / Cp PT — vs Fortran (scalars only)
    # =========================================================================
    cv_pt = ar_model.cv_residual_pt(n, p, t, root="stable")
    assert np.isclose(cv_pt, cv_pt_f), f"Cv_pt {name}"

    cp_pt = ar_model.cp_residual_pt(n, p, t, root="stable")
    assert np.isclose(cp_pt, cp_pt_f), f"Cp_pt {name}"

    # =========================================================================
    # Ge — vs Fortran + individual derivatives
    # =========================================================================
    ge, ge_d = ar_model.gibbs_excess(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.isclose(ge, ge_f), f"Ge {name}"
    assert np.isclose(ge_d["dt"], get_f), f"GeT {name}"
    assert np.isclose(ge_d["dp"], gep_f), f"GeP {name}"
    assert np.allclose(ge_d["dn"], gen_f), f"Gen {name}"

    # Individual
    ge_i = ar_model.gibbs_excess(n, p, t, root="stable")
    assert np.isclose(ge_i, ge_f), f"Ge individual {name}"

    _, d = ar_model.gibbs_excess(n, p, t, root="stable", dt=True)
    assert np.isclose(d["dt"], get_f), f"GeT individual {name}"

    _, d = ar_model.gibbs_excess(n, p, t, root="stable", dp=True)
    assert np.isclose(d["dp"], gep_f), f"GeP individual {name}"

    _, d = ar_model.gibbs_excess(n, p, t, root="stable", dn=True)
    assert np.allclose(d["dn"], gen_f), f"Gen individual {name}"

    # =========================================================================
    # Ve — vs Fortran (scalar only)
    # =========================================================================
    ve = ar_model.volume_excess(n, p, t, root="stable")
    assert np.isclose(ve, ve_f), f"Ve {name}"

    # =========================================================================
    # lngamma — vs Fortran + individual derivatives
    # =========================================================================
    lngamma, lngamma_d = ar_model.ln_gamma(
        n, p, t, root="stable", dt=True, dp=True, dn=True
    )

    assert np.allclose(lngamma, lngamma_f), f"lngamma {name}"
    assert np.allclose(lngamma_d["dt"], dlngammadt_f), f"dlngammadT {name}"
    assert np.allclose(lngamma_d["dp"], dlngammadp_f), f"dlngammadP {name}"
    assert np.allclose(lngamma_d["dn"], dlngammadn_f), f"dlngammadn {name}"

    # Individual
    lngamma_i = ar_model.ln_gamma(n, p, t, root="stable")
    assert np.allclose(lngamma_i, lngamma_f), f"lngamma individual {name}"

    _, d = ar_model.ln_gamma(n, p, t, root="stable", dt=True)
    assert np.allclose(d["dt"], dlngammadt_f), f"dlngammadT individual {name}"

    _, d = ar_model.ln_gamma(n, p, t, root="stable", dp=True)
    assert np.allclose(d["dp"], dlngammadp_f), f"dlngammadP individual {name}"

    _, d = ar_model.ln_gamma(n, p, t, root="stable", dn=True)
    assert np.allclose(d["dn"], dlngammadn_f), f"dlngammadn individual {name}"

    # =========================================================================
    # Ae, He, Ue, Se — vs Fortran (scalars only)
    # =========================================================================
    ae = ar_model.helmoltz_excess(n, p, t, root="stable")
    assert np.isclose(ae, ae_f), f"Ae {name}"

    he = ar_model.enthalpy_excess(n, p, t, root="stable")
    assert np.isclose(he, he_f), f"He {name}"

    ue = ar_model.internal_energy_excess(n, p, t, root="stable")
    assert np.isclose(ue, ue_f), f"Ue {name}"

    se = ar_model.entropy_excess(n, p, t, root="stable")
    assert np.isclose(se, se_f), f"Se {name}"
