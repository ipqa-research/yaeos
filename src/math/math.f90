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

   function derivative_dxk_dni(n) result(dxk_dni)
      !! # derivative_dxk_dni
      !!
      !! # Description
      !! Calculate the mole fraction first derivatives respect to mole numbers
      real(pr), intent(in) :: n(:)
      real(pr) :: dxk_dni(size(n), size(n))

      real(pr) :: n_tot
      integer :: nc, k, i

      n_tot = sum(n)
      nc = size(n)

      dxk_dni = 0.0_pr
      do concurrent(k=1:nc, i=1:nc)
         if (k == i) then
            dxk_dni(k,i) = (n_tot - n(i)) / n_tot**2
         else
            dxk_dni(k,i) = -n(k) / n_tot**2
         end if
      end do
   end function derivative_dxk_dni

   function derivative_d2xk_dnidnj(n) result(d2xk_dnidnj)
      !! # derivative_d2xk_dnidnj
      !!
      !! # Description
      !! Calculate the mole fraction second derivatives respect to mole numbers
      real(pr), intent(in) :: n(:)
      real(pr) :: d2xk_dnidnj(size(n), size(n), size(n))

      real(pr) :: n_tot
      integer :: nc, k, i, j

      n_tot = sum(n)
      nc = size(n)

      d2xk_dnidnj = 0.0_pr
      do concurrent (k=1:nc, i=1:nc, j=1:nc)
         if (i==k .and. j==k) then
            d2xk_dnidnj(k,i,j) = -2 * (n_tot - n(i)) / n_tot**3
         else if (i==k) then
            d2xk_dnidnj(k,i,j) = (2 * n(i) - n_tot) / n_tot**3
         else if (j==k) then
            d2xk_dnidnj(k,i,j) = (2 * n(j) - n_tot) / n_tot**3
         else
            d2xk_dnidnj(k,i,j) = 2 * n(k) / n_tot**3
         end if
      end do
   end function derivative_d2xk_dnidnj

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
