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
   use yaeos__math_linalg, only: solve_system

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

      integer :: i, j, nc

      nc = size(x)
      dn = 0

      forall(i=1:nc, j=1:nc, i /= j)
         dn(i) = dn(i) - x(j) * dx(j)
      end forall

      dn = dx - dn
   end function
end module yaeos__math
