module yaeos_substance
    !! Ya_EoS Subtance module.
    !!
    !! Module containing pure components properties and parameters.
    use yaeos_constants, only: pr

    type :: Substances
        !! Set of pure components
        character(len=50), allocatable :: names(:) !! Composition names.
        real(pr), allocatable :: tc(:) !! Critical Temperature [K]
        real(pr), allocatable :: pc(:) !! Critical Pressure [bar]
        real(pr), allocatable :: w(:) !! Acentric factor
    end type
end module