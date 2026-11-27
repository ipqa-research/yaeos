module yaeos__substance
    !! # Component properties
    !! Derived type for pure-component data shared by thermodynamic models.
    !!
    !! # Description
    !! `Substances` stores component names and critical properties in the order
    !! used by a model's composition vectors. Optional allocatable arrays may
    !! be left unallocated when a model does not require those properties.
    !!
    !! # Examples
    !! ```fortran
    !! type(Substances) :: components
    !! components = Substances(tc=Tc, pc=Pc, w=omega)
    !! ```
    use yaeos__constants, only: pr

    type :: Substances
        !! # Substances
        !! Pure-component identifiers and critical properties.
        !!
        !! # Description
        !! Arrays use the model's component ordering and are allocated by
        !! assignment when initialized from constructor arguments.
        !!
        !! # Examples
        !! ```fortran
        !! components = Substances(names=names, tc=Tc, pc=Pc)
        !! ```
        character(len=50), allocatable :: names(:) !! Composition names.
        real(pr), allocatable :: tc(:) !! Critical Temperature [K]
        real(pr), allocatable :: pc(:) !! Critical Pressure [bar]
        real(pr), allocatable :: w(:) !! Acentric factor
        real(pr), allocatable :: vc(:) !! Critical Volume [L/mol]
    end type
end module
