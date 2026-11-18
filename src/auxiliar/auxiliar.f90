module yaeos__auxiliar
   !! # General auxiliary utilities
   !! Small convenience procedures shared by model, solver, and phase-equilibrium code.
   !!
   !! # Description
   !! The generic `optval` returns either a present optional scalar or a supplied
   !! default. `sort` sorts a real array in place using `stdlib_sorting`.
   !!
   !! # Examples
   !! ```fortran
   !! niter = optval(max_iters, 100)
   !! call sort(values)
   !! ```
   use yaeos__constants, only: pr
   implicit none

   interface optval
      module procedure optval_integer, optval_real, optval_character, optval_logical
   end interface optval

contains

   integer function optval_integer(val, default)
      !! # optval_integer
      !! Return the optional integer value or its fallback.
      !!
      !! # Description
      !! Selects `val` when present and otherwise returns `default`.
      !!
      !! # Examples
      !! ```fortran
      !! iterations = optval(max_iters, 100)
      !! ```
      integer, optional, intent(in) :: val !! Optional integer value.
      integer, intent(in) :: default !! Value returned when `val` is absent.

      if (present(val)) then
         optval_integer = val
      else
         optval_integer = default
      end if
   end function optval_integer
   
   logical function optval_logical(val, default)
      !! # optval_logical
      !! Return the optional logical value or its fallback.
      !!
      !! # Description
      !! Selects `val` when present and otherwise returns `default`.
      logical, optional, intent(in) :: val !! Optional logical value.
      logical, intent(in) :: default !! Value returned when `val` is absent.

      if (present(val)) then
         optval_logical = val
      else
         optval_logical = default
      end if
   end function optval_logical

   real(pr) function optval_real(val, default)
      !! # optval_real
      !! Return the optional real value or its fallback.
      !!
      !! # Description
      !! Selects `val` when present and otherwise returns `default`.
      real(pr), optional, intent(in) :: val !! Optional real value.
      real(pr), intent(in) :: default !! Value returned when `val` is absent.

      if (present(val)) then
         optval_real = val
      else
         optval_real = default
      end if
   end function optval_real

   function optval_character(val, default)
      !! # optval_character
      !! Return the optional character value or its fallback.
      !!
      !! # Description
      !! Selects `val` when present and otherwise returns `default`; the result
      !! is allocated to the selected string length.
      character(len=*), optional, intent(in) :: val !! Optional character value.
      character(len=*), intent(in) :: default !! Value returned when `val` is absent.
      character(len=:), allocatable :: optval_character

      if (present(val)) then
         optval_character = val
      else
         optval_character = default
      end if
   end function optval_character

   subroutine sort(array, idx)
      !! # sort
      !! Sort a real array in ascending order in place.
      !!
      !! # Description
      !! Delegates the sorting operation to `stdlib_sorting`. The optional `idx`
      !! argument is present in the interface but is not currently populated.
      !!
      !! # Examples
      !! ```fortran
      !! call sort(values)
      !! ```
      use stdlib_sorting, only: std => sort
      real(pr), intent(in out) :: array(:) !! Values to sort; overwritten in ascending order.
      integer, optional, intent(out) :: idx(:) !! Reserved optional index output; currently unchanged.

      call std(array)

   end subroutine sort
end module yaeos__auxiliar
