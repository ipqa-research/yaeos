module constants
   use iso_fortran_env, only: real32, real64, real128
   implicit none

   integer, parameter :: pr = real64
   real(pr), parameter :: R = 0.08314472
   character(len=254) :: database_path = "database/"
   character(len=1) :: path_sep = "/"
end module constants