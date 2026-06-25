program main
   use yaeos, only: pr
   use yaeos__newton_solver, only: newton
   use testing_aux, only: assert, test_title
   use auxiliar_functions, only: allclose

   real(pr) :: x(2)
   real(pr) :: Fval(2), dFval(2, 2)

   write(*, *) test_title("LLM NEWTON")

   x = [1., 39.]
   call newton(F, x)
   call F(X, Fval, dFval)

   call assert(norm2(Fval) < 1e-9_pr, "Should solve the equation")
   call assert(allclose(X , [1._pr, 3._pr], 1e-9_pr), "Correct solution")

contains

   subroutine F(X, Fval, dF)
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: Fval(:)
      real(pr), intent(out) :: dF(:, :)

      Fval(1) = X(1)**2 + X(2)**2 - 10
      Fval(2) = X(1) - 1

      dF(1, 1) = 2*X(1)
      dF(1, 2) = 2*X(2)
      dF(2, 1) = 1
   end subroutine F
end program main
