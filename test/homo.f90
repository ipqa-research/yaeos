program main
   use yaeos__math_nonlinearsolvers, only: homotopy, newton
   implicit none
   real(8) :: X(3), F(3), J(3,3)
   integer :: its

   X = [0.5_8, 0.5_8, 0.5_8]


   ! call newton(sub=foo, x=X, tol=1e-8_8, max_its=50, its=its)
   call homotopy(sub=foo, x=X, tol=1e-8_8, max_its=50, its=its)

   call foo(X, F, J)
   print *, "Final F:", F

contains

   subroutine foo(X, F, J)
      real(8), intent(in) :: X(:)
      real(8), intent(out) :: F(:)
      real(8), intent(out) :: J(:, :)

      F(1) = X(1)**2 + X(2)**2 + X(3)**2 - 1.0_8
      F(2) = X(1) + X(2) - X(3)
      F(3) = X(1) - X(2) + X(3)

      J(1,1) = 2.0_8*X(1)
      J(1,2) = 2.0_8*X(2)
      J(1,3) = 2.0_8*X(3)

      J(2,1) = 1.0_8
      J(2,2) = 1.0_8
      J(2,3) = -1.0_8

      J(3,1) = 1.0_8
      J(3,2) = -1.0_8
      J(3,3) = 1.0_8

   end subroutine
end program