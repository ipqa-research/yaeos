module yaeos_constants
   !! Constants used on the whole package
   use iso_fortran_env, only: real32, real64, real128
   
   implicit none

   integer, parameter :: pr = real64 !! Used precision
   real(pr), parameter :: R = 0.08314462618_pr !! Ideal Gas constant
   character(len=254) :: database_path = "database" !! Path to find database
   character(len=1) :: path_sep = "/" !! File separator (to preprocess on Win or Mac/linux)

   real(pr), parameter :: NOT_IMPLEMENTED = huge(R)
end module