module yaeos__models_base
    !! Basic element of a thermodynamic model.
    use yaeos__constants, only: pr
    use yaeos__substance, only: Substances
    implicit none

    type, abstract :: BaseModel
        !! Base model type.
        !!
        !! Contains the important parts of most models and other procedures.
        type(Substances) :: components 
            !! Substances contained in the module
    end type

    type :: PropnVT(nc)
        !! Thermodynamic property.
        integer, len :: nc
        real(pr) :: val !! Value.
        real(pr) :: dv !! First derivative with volume.
        real(pr) :: dv2 !! Second derivative with volume.
        real(pr) :: dt !! First derivative with temperature.
        real(pr) :: dt2 !! Second derivative with temperature.
        real(pr) :: dtv !! Crossed derivative with temperature and volume.
        real(pr) :: dn(nc) !! Derivative with mole number.
        real(pr) :: dvn(nc) !! Crossed derivative with volume and moles number.
        real(pr) :: dtn(nc) !! Crossed derivative with temperature and moles number.
        real(pr) :: dn2(nc, nc) !! Second derivative with mole number.
    end type
end module