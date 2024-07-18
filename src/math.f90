module yaeos__math
   !! # Mathematical methods for `yaeos`
   !!
   !! # Description
   !! This module provides all the relevant mathematical functions used in this
   !! library. Most important ones are:
   !!
   !! - newton: Newton solving method
   !! - solve_system: Solving linear system Ax = b
   !! - continuation: Continuation method for line tracing
   !!
   !! # Examples
   !!
   !! ## Squared error calculation
   !! ```fortran
   !!  use yaeos__math, only: sq_error
   !!  real(pr) :: x = 2.5, y = 3.0, error
   !!  print *, sq_error(2.5, 3.0)
   !! ------------------------------------
   !! ```
   !!
   !! ```fortran
   !!  use yaeos__math, only: sq_error
   !!  real(pr) :: x = [2.5, 5.0], y = [3.0, 4.5], error
   !!  ! It also works with arrays
   !!  print *, sq_error(x, y)
   !! ```

   use yaeos__math_continuation, only: continuation
   use yaeos__math_linalg, only: solve_system, cubic_roots
   use yaeos__constants, only: pr

   implicit none

   abstract interface
      subroutine f_1d(x, f, df)
         import pr
         real(pr), intent(in) :: x
         real(pr), intent(out) :: f
         real(pr), intent(out) :: df
      end subroutine f_1d
   end interface

   interface newton
      module procedure :: newton_1d
   end interface newton

contains
   elemental real(pr) function sq_error(exp, pred)
      !! # Squared error between two values.
      !!
      !! # Description
      !! ...
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  error = sq_error(true_value, model_value)
      !! ```
      use yaeos__constants, only: pr
      real(pr), intent(in) :: exp
      real(pr), intent(in) :: pred
      sq_error = ((exp - pred)/exp)**2
   end function sq_error

   function dx_to_dn(x, dx) result(dn)
      !! # dx_to_dn
      !!
      !! # Description
      !! Convert the mole fraction derivatives of a quantity (calculated
      !! so they do not sum to 1) to mole number derivatives (where the mole
      !! fractions do sum to one). Requires the derivatives and the mole fractions
      !! of the mixture.
      !! From [https://chemicals.readthedocs.io/chemicals.utils.html?highlight=dxs_to_dns#chemicals.utils.dxs_to_dns](Chemicals (Python))
      use yaeos__constants, only: pr
      real(pr), intent(in) :: x(:)
      real(pr), intent(in) :: dx(:)
      real(pr) :: dn(size(x))

      real(pr) :: sum_xdx

      dn = 0

      sum_xdx = sum(x * dx)

      dn = dx - sum_xdx
   end function dx_to_dn

   subroutine newton_1d(f, x, tol, max_iters)
      procedure(f_1d) :: f
      real(pr), intent(in out) :: x
      real(pr), intent(in) :: tol
      integer, intent(in) :: max_iters

      integer :: i
      real(pr) :: fval, df, step


      fval = 10
      step = 10

      do i=1, max_iters
         if (abs(fval) < tol .or. abs(step) < tol)  exit
         call f(x, fval, df)

         step = fval/df

         do while (abs(step) > 0.5 * abs(x))
            step = step/2
         end do

         x = x - step
      end do
   end subroutine newton_1d
end module yaeos__math
