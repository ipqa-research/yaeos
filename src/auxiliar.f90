module yaeos__auxiliar
   use yaeos__constants, only: pr
   implicit none

   interface optval
      module procedure optval_integer, optval_real
   end interface optval

contains

   integer function optval_integer(val, default)
      !! Set a value to a default if it is not defined
      use stdlib_optval, only: std => optval
      integer, optional, intent(in out) :: val
      integer, intent(in) :: default
      optval_integer = std(val, default)
   end function optval_integer

   real(pr) function optval_real(val, default)
      !! Set a value to a default if it is not defined
      use stdlib_optval, only: std => optval
      real(pr), optional, intent(in out) :: val
      real(pr), intent(in) :: default

      optval_real = std(val, default)
   end function optval_real

   subroutine sort(array, idx)
      use stdlib_sorting, only: std => sort
      !! Sort an array and return the indexes
      real(pr), intent(in out) :: array(:)
      integer, optional, intent(out) :: idx(:)

      call std(array)

   end subroutine sort
end module yaeos__auxiliar
