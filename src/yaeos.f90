module yaeos
   !! # yaeos thermodynamic library
   !! Public top-level module for EOS properties and phase-equilibrium calculations.
   !!
   !! # Description
   !! This module re-exports the shared constants, component data, model types,
   !! consistency tools, and phase-equilibrium procedures. Most users can start
   !! with `use yaeos` and then construct a model, evaluate properties, and run
   !! equilibrium routines without importing implementation submodules.
   !!
   !! The main API groups are:
   !!
   !! - [[yaeos__constants(module)]]: Units, real precision, and shared constants.
   !! - [[yaeos__consistency(module)]]: Numerical consistency checks for models.
   !! - [[yaeos__substance(module)]]: Component names and pure-fluid properties.
   !! - [[yaeos__models(module)]]: EOS and excess-Gibbs models and base types.
   !! - [[yaeos__equilibria(module)]]: Stability, saturation, flash, and critical routines.
   !!
   !! # Examples
   !! ```fortran
   !! use yaeos
   !! class(ArModel), allocatable :: model
   !! type(EquilibriumState) :: result
   !! real(pr) :: n(2)
   !! real(pr) :: z(2)
   !! real(pr) :: Tc(2)
   !! real(pr) :: Pc(2)
   !! real(pr) :: omega(2)
   !! real(pr) :: V
   !! real(pr) :: T
   !! real(pr) :: P
   !! model = PengRobinson76(Tc, Pc, omega)
   !! call model%pressure(n, V, T, P)
   !! z = n / sum(n)
   !! result = flash(model, z, T, p_spec=P)
   !! ```
   use yaeos__constants
   use yaeos__consistency
   use yaeos__substance
   use yaeos__models
   use yaeos__equilibria
   character(len=*), parameter :: version="4.5.4" !! This version.
end module yaeos
