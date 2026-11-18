module yaeos__constants
   !! # Library constants and shared settings
   !! Precision, gas constant, unit conventions, volume-root identifiers, and
   !! shared runtime settings used throughout `yaeos`.
   !!
   !! # Description
   !! `pr` selects the library real kind and `R` defines the default unit
   !! system (bar, litre, kelvin, mole). `root_kinds` names the stable, liquid,
   !! and vapor roots used by EOS volume calculations.
   !!
   !! # Examples
   !! ```fortran
   !! use yaeos, only: pr, R, root_kinds
   !! real(pr) :: temperature = 300.0_pr
   !! print *, R, root_kinds%stable
   !! ```
   use iso_fortran_env, only: real32, real64, real128
   
   implicit none

   integer, parameter :: pr = real64 !! Used precision
   real(pr), parameter :: R = 0.08314462618_pr !! Ideal Gas constant [bar L / (mol K)]
   character(len=254) :: database_path = "database" !! Path to find database
   character(len=1) :: path_sep = "/" !! File separator (to preprocess on Win or Mac/linux)

   real(pr), parameter :: NOT_IMPLEMENTED = huge(R) !! Sentinel used for unavailable property values.
   logical :: solving_volume = .false. !! Internal guard indicating an active volume solve.
   
   type :: KindEnum
      !! # KindEnum
      !! Named root selectors used by EOS volume calculations.
      !!
      !! # Description
      !! The integer fields identify the stable, liquid, and vapor roots for
      !! routines accepting a `root_type` selector.
      !!
      !! # Examples
      !! ```fortran
      !! call model%lnphi_pt(n, P, T, root_type="liquid", lnPhi=lnphi)
      !! ```
      integer :: stable=0 !! Select the thermodynamically stable volume root.
      integer :: liquid=1 !! Select the liquid volume root.
      integer :: vapor=2 !! Select the vapor volume root.
   end type KindEnum

   type(KindEnum), parameter :: root_kinds = KindEnum() !! KindEnum instance

end module
