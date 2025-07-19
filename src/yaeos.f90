module yaeos
   !! Yet Another Equation-Of-State (library)
   !!
   !! Library to use EoS-based calculations. This main module imports all the
   !! relevant constants, procedures and objects to have better access to them
   !! The main submodules that it uses are:
   !!
   !! - [[yaeos__constants(module)]]: All the relevant costants and also the used precision (default=double precision).
   !! - [[yaeos__consistency(module)]]: Tools to evalaute the consistency of Ar and Ge models.
   !! - [[yaeos__substance(module)]]: Derived type that holds the important data (for example, critical constants) from a mixture.
   !! - [[yaeos__models(module)]]: All the implemented models, also their base types for making extensions.
   !! - [[yaeos__equilibria(module)]]: Phase equilibria related procedures.
   use yaeos__constants
   use yaeos__consistency
   use yaeos__substance
   use yaeos__models
   use yaeos__equilibria
   character(len=*), parameter :: version="4.2.2" !! This version.
end module
