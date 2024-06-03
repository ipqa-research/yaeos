module yaeos__math
   !! Mathematical methods for `yaeos`
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
   !! ```fortran
   !!  A basic code example
   !! ```

   use yaeos__math_continuation, only: continuation
   use yaeos__math_linalg, only: solve_system

contains
   real(pr) function sq_error(exp, pred)
      !! Squared error between two values.
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
end module yaeos__math
