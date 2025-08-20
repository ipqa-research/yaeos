module yaeos__auxiliar
   !! # Auxiliary Functions and Utilities
   !!
   !! This module provides utility functions used throughout yaeos,
   !! particularly for handling optional arguments and providing
   !! default values.
   !!
   !! ## Main Features
   !!
   !! - **Optional value handling**: `optval` interface for setting defaults
   !! - **Type-generic**: Works with integers, reals, and characters
   !!
   !! ## Usage Examples
   !!
   !! ```fortran
   !! use yaeos__auxiliar, only: optval
   !!
   !! subroutine my_procedure(required_arg, optional_arg)
   !!     real(pr), intent(in) :: required_arg
   !!     real(pr), optional, intent(in) :: optional_arg
   !!     
   !!     real(pr) :: tol
   !!     
   !!     ! Set default tolerance if not provided
   !!     tol = optval(optional_arg, 1e-6_pr)
   !! end subroutine
   !! ```
   !!
   !! The `optval` interface provides a clean way to handle optional
   !! arguments with default values, improving code readability and
   !! reducing boilerplate code.
   use yaeos__constants, only: pr
   implicit none

   interface optval
      module procedure optval_integer, optval_real, optval_character
   end interface optval

contains

   integer function optval_integer(val, default)
      !! Set a value to a default if it is not defined
      integer, optional, intent(in) :: val
      integer, intent(in) :: default

      if (present(val)) then
         optval_integer = val
      else
         optval_integer = default
      end if
   end function optval_integer

   real(pr) function optval_real(val, default)
      !! Set a value to a default if it is not defined
      real(pr), optional, intent(in) :: val
      real(pr), intent(in) :: default

      if (present(val)) then
         optval_real = val
      else
         optval_real = default
      end if
   end function optval_real

   function optval_character(val, default)
      !! Set a value to a default if it is not defined
      character(len=*), optional, intent(in) :: val
      character(len=*), intent(in) :: default
      character(len=:), allocatable :: optval_character

      if (present(val)) then
         optval_character = val
      else
         optval_character = default
      end if
   end function optval_character

   subroutine sort(array, idx)
      use stdlib_sorting, only: std => sort
      !! Sort an array and return the indexes
      real(pr), intent(in out) :: array(:)
      integer, optional, intent(out) :: idx(:)

      call std(array)

   end subroutine sort
end module yaeos__auxiliar
