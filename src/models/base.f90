module yaeos__models_base
    !! # Base thermodynamic model
    !! Shared base type for all thermodynamic model families.
    !!
    !! # Description
    !! `BaseModel` stores component properties through its `components` field.
    !! Residual Helmholtz models and excess-Gibbs models extend this base type
    !! and add their respective calculation interfaces.
    !!
    !! # Examples
    !! ```fortran
    !! class(BaseModel), allocatable :: model
    !! model = PengRobinson76(Tc, Pc, omega)
    !! ```
    use yaeos__substance, only: Substances
    implicit none

    type, abstract :: BaseModel
        !! # BaseModel
        !! Common parent type for models that operate on a component set.
        !!
        !! # Description
        !! Provides component metadata used by property and phase-equilibrium
        !! routines. Concrete model families extend this type with their
        !! thermodynamic equations and model-specific parameters.
        !!
        !! # Examples
        !! See the `ArModel` and `GeModel` interfaces for concrete extensions.
        type(Substances) :: components
        !! Pure-component properties for the model's component order.
    end type
end module
