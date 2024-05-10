module yaeos__math_continuation
   use yaeos_constants, only: pr
   use yaeos__math_linalg, only: solve_system

   abstract interface
      subroutine continuation_function(X, ns, S, F, dF, dFdS)
         import pr
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)
      end subroutine continuation_function

      subroutine process(X, ns, S, dS, dXdS, iterations)
         import pr
         real(pr), intent(in out) :: X(:)
         integer, intent(in out) :: ns
         real(pr), intent(in out) :: S
         real(pr), intent(in out) :: dS
         real(pr), intent(in out) ::  dXdS(:)
         integer, intent(in) :: iterations
      end subroutine process

      subroutine continuation_solver(fun, iters, X, ns, S, max_iters, F, dF, dFdS, tol &
         )
         import pr
         procedure(foo) :: fun !! Function to solve
         integer,  intent(out)    :: iters !! Number of iterations needed
         real(pr), intent(in out) :: X(:)  !! Variables vector
         integer, intent(in) :: ns !! Specification number
         real(pr), intent(in) :: S !! Specification value
         integer, intent(in)      :: max_iters !! Maximum iterations
         real(pr), intent(out)    :: F(:) !! Function values at solved point
         real(pr), intent(out)    :: df(:, :) !! Jacobian values
         real(pr), intent(out)    :: dfds(:) !! dFdS
         real(pr), intent(in) :: tol !! Solver tolerance
      end subroutine continuation_solver
   end interface

   abstract interface
      subroutine foo(X, ns, S, F, dF, dFdS)
         import pr
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)
      end subroutine foo
   end interface

contains

   function continuation(&
      f, X0, ns0, S0, dS0, max_points, solver_tol, &
      update_specification, postprocess, solver &
      ) result(XS)
      procedure(continuation_function) :: f
      real(pr), intent(in) :: X0(:)
      integer, intent(in) :: ns0
      real(pr), intent(in) :: S0
      real(pr), intent(in) :: dS0
      integer, intent(in) :: max_points
      real(pr), intent(in) :: solver_tol
      procedure(process), optional :: update_specification
      procedure(process), optional :: postprocess
      procedure(continuation_solver), optional :: solver
      real(pr) :: XS(max_points, size(X0))

      real(pr) :: X(size(X0)), S, fval(size(X0)), dF(size(X0), size(X0)), dFdS(size(X0))
      real(pr) :: dXdS(size(X0))
      integer :: ns
      real(pr) :: dS

      integer :: i, newton_its

      integer :: max_iters = 100

      X = X0
      ns = ns0
      dS = dS0
      S = S0

      XS = 0

      do i=1,max_points

         if (present(solver)) then
            call solver(f, newton_its, X, ns, S, max_iters, fval, dF, dFdS, solver_tol)
         else
            call full_newton(&
               f, newton_its, X, ns, S, max_iters, fval, dF, dFdS, solver_tol &
               )
         end if
         if (newton_its >= max_iters) exit

         XS(i, :) = X

         dXdS = solve_system(dF, -dFdS)

         if (present(update_specification)) then
            call update_specification(X, ns, S, dS, dXdS, newton_its)
         else
            ns = maxloc(abs(dXdS), dim=1)
            dS = dXdS(ns)*dS
            dXdS = dXdS/dXdS(ns)
            dS = sign(minval(abs([0.05_pr, dS])), dS)
         end if

         if (present(postprocess)) then
            call postprocess(X, ns, S, dS, dXdS, newton_its)
         end if

         X = X + dXdS * dS
         S = X(ns)
      end do
   end function continuation

   subroutine full_newton(&
      fun, iters, X, ns, S, max_iters, F, dF, dFdS, tol &
      )
      !! Subroutine to solve a point.
      !!
      !! Procedure that solves a point with the Newton-Raphson method.
      use stdlib_optval, only: optval
      use yaeos__math_linalg, only: solve_system
      procedure(foo) :: fun !! Function to solve
      integer,  intent(out)    :: iters !! Number of iterations needed
      real(pr), intent(in out) :: X(:)  !! Variables vector
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      integer, intent(in)      :: max_iters !! Maximum iterations
      real(pr), intent(out)    :: F(:) !! Function values at solved point
      real(pr), intent(out)    :: df(:, :) !! Jacobian values
      real(pr), intent(out)    :: dfds(:) !! dFdS
      real(pr), optional, intent(in) :: tol

      real(pr) :: dX(size(X)), solve_tol

      solve_tol = optval(tol, 1.e-5_pr)

      dX = 20
      F = 500
      newton: do iters = 1, max_iters
         ! Converged point
         if (maxval(abs(dx/x)) < solve_tol .or. maxval(abs(F)) < solve_tol) exit newton

         call fun(X, ns, S, F, dF, dFdS)

         dX = solve_system(dF, -F)

         do while(maxval(abs(dx)) > 1)
            dX = dX/2
         end do

         X = X + dX
      end do newton
   end subroutine full_newton
end module yaeos__math_continuation
