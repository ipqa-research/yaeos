module yaeos__constants
   !! # yaeos Constants and Global Settings
   !!
   !! This module defines all the fundamental constants used throughout
   !! the yaeos library, as well as global configuration parameters.
   !!
   !! ## Precision Control
   !!
   !! The working precision for all real numbers is controlled by the
   !! `pr` parameter, which defaults to double precision (real64).
   !!
   !! ```fortran
   !! real(pr) :: my_variable  ! Will be double precision
   !! ```
   !!
   !! ## Physical Constants
   !!
   !! The ideal gas constant R is defined with units consistent with
   !! the yaeos unit system:
   !!
   !! ```fortran
   !! R = 0.08314462618 bar⋅L/(mol⋅K)
   !! ```
   !!
   !! ## Unit System
   !!
   !! yaeos uses the following consistent unit system:
   !! - Temperature: Kelvin [K]
   !! - Pressure: bar [bar]
   !! - Volume: Liters [L]
   !! - Amount: moles [mol]
   !! - Energy: bar⋅L (= 0.1 J)
   !!
   !! ## Configuration Parameters
   !!
   !! - `database_path`: Path to parameter databases
   !! - `path_sep`: File path separator (OS-dependent)
   !!
   !! ## Enumeration Types
   !!
   !! The `KindEnum` type provides standardized phase identifiers:
   !! - `stable`: Most stable phase (0)
   !! - `liquid`: Liquid phase (1)  
   !! - `vapor`: Vapor phase (2)
   !!
   !! ## Usage Examples
   !!
   !! ```fortran
   !! use yaeos__constants
   !!
   !! ! Use working precision
   !! real(pr) :: temperature = 298.15_pr
   !!
   !! ! Use gas constant
   !! real(pr) :: ideal_pressure = n * R * T / V
   !!
   !! ! Use phase identifiers
   !! integer :: phase_type = root_kinds%vapor
   !! ```
   use iso_fortran_env, only: real32, real64, real128
   
   implicit none

   integer, parameter :: pr = real64 !! Used precision
   real(pr), parameter :: R = 0.08314462618_pr !! Ideal Gas constant [bar L / (mol K)]
   character(len=254) :: database_path = "database" !! Path to find database
   character(len=1) :: path_sep = "/" !! File separator (to preprocess on Win or Mac/linux)

   real(pr), parameter :: NOT_IMPLEMENTED = huge(R)
   logical :: solving_volume = .false.
   
   type :: KindEnum
      !! Enumeration of the possible phases that can be used
      !! in the envelope calculations
      integer :: stable=0
      integer :: liquid=1
      integer :: vapor=2
   end type KindEnum

   type(KindEnum), parameter :: root_kinds = KindEnum() !! KindEnum instance

end module