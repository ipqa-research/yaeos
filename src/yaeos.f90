module yaeos
   !! Yet Another Equation-Of-State (library)
   !!
   !! Library to use EoS-based calculations. This main module imports all the
   !! relevant constants, procedures and objects to have better access to them
   !! The main submodules that it uses are:
   !!
   !! - [[yaeos_constants(module)]]: All the relevant costants and also the used precision (default=double precision).
   !! - [[yaeos_consistency(module)]]: Tools to evalaute the consistency of Ar and Ge models.
   !! - [[yaeos_substance(module)]]: Derived type that holds the important data (for example, critical constants) from a mixture.
   !! - [[yaeos_models(module)]]: All the implemented models, also their base types for making extensions.
   !! - [[yaeos_thermoprops(module)]]: Available thermodynamic properties to calculate.
   !! - [[yaeos_equilibria(module)]]: Phase equilibria related procedures.
   use yaeos_constants
   use yaeos_consistency
   use yaeos_substance
   use yaeos_models
   use yaeos_thermoprops
   use yaeos_equilibria
   character(len=*), parameter :: version="0.2.0b1" !! This version.
end module
