module yaeos__constants
   !! Constants used on the whole package
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