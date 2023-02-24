module models
    use ar_models

    private

    public :: ArModel, residual_helmholtz
    public :: size, alloc

    public :: setup_model

    ! Cubic EoS
    public :: attractive_parameter
    public :: repulsive_parameter
    public :: del1_parameter
    public :: del2_parameter
    public :: ar

    public :: CubicEOS
    public :: PR76

    ! Mixing rules
    public :: CubicMixingRule, BinaryMixingRule
    
    type :: Model
        class(ArModel), pointer :: residual
    end type Model
end module models