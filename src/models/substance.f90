module yaeos__substance
    !! # Pure Component Properties Module
    !!
    !! This module defines the data structures for storing pure component
    !! properties and parameters used in thermodynamic models.
    !!
    !! ## Main Features
    !!
    !! - **Component database**: Storage for pure component properties
    !! - **Critical constants**: Tc, Pc, Vc, acentric factor
    !! - **Component identification**: Names and identifiers
    !! - **Extensible design**: Easy to add new properties
    !!
    !! ## Data Structure
    !!
    !! The `Substances` type contains arrays of component properties:
    !!
    !! ```fortran
    !! type(Substances) :: components
    !! 
    !! ! Properties available:
    !! components%names(:)  ! Component names
    !! components%tc(:)     ! Critical temperatures [K]
    !! components%pc(:)     ! Critical pressures [bar]  
    !! components%w(:)      ! Acentric factors [-]
    !! components%vc(:)     ! Critical volumes [L/mol]
    !! ```
    !!
    !! ## Usage Examples
    !!
    !! ### Manual Component Setup
    !! ```fortran
    !! use yaeos__substance
    !! 
    !! type(Substances) :: comps
    !! 
    !! ! Allocate for 2 components
    !! allocate(comps%names(2), comps%tc(2), comps%pc(2), comps%w(2))
    !! 
    !! ! Ethane + n-Butane system
    !! comps%names = ["ethane   ", "n-butane "]
    !! comps%tc = [305.32, 425.12]  ! K
    !! comps%pc = [48.72, 37.96]    ! bar
    !! comps%w = [0.0995, 0.2002]   ! dimensionless
    !! ```
    !!
    !! ### Integration with Models
    !! ```fortran
    !! class(ArModel), allocatable :: model
    !! 
    !! ! Most models accept component properties directly
    !! model = PengRobinson76(comps%tc, comps%pc, comps%w)
    !! 
    !! ! Or can be stored within model structures
    !! model%components = comps
    !! ```
    !!
    !! ## Database Integration
    !!
    !! This type is designed to integrate with component databases
    !! and parameter files for convenient setup of complex systems.
    !!
    !! ## See Also
    !!
    !! - [[yaeos__models_base(module)]] - Base model containing substances
    !! - Component databases in the `database/` directory
    use yaeos__constants, only: pr

    type :: Substances
        !! Set of pure components
        character(len=50), allocatable :: names(:) !! Composition names.
        real(pr), allocatable :: tc(:) !! Critical Temperature [K]
        real(pr), allocatable :: pc(:) !! Critical Pressure [bar]
        real(pr), allocatable :: w(:) !! Acentric factor
        real(pr), allocatable :: vc(:) !! Critical Volume [L/mol]
    end type
end module