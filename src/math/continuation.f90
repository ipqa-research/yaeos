module yaeos__math_continuation
   !! Implementation of Algower's numerical continuation method.
   use yaeos_constants, only: pr
   use yaeos__math_linalg, only: solve_system

   implicit none

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
   end interface

   abstract interface
      subroutine process(X, ns, S, dS, dXdS, iterations)
         !! Subroutine to make variation in the method after a point converged
         import pr
         real(pr), intent(in out) :: X(:) !! Vector of variables \(X\)
         integer, intent(in out) :: ns !! Position of specified variable
         real(pr), intent(in out) :: S !! Specification variable value
         real(pr), intent(in out) :: dS !! Step of specification in the method
         real(pr), intent(in out) ::  dXdS(:) !! \(\frac{dX}{dS}\)
         integer, intent(in) :: iterations !! Iterations needed to converge point
      end subroutine process

      logical function continuation_stopper(X, ns, S, dS, dXdS, iterations)
         !! Function that returns true if the method should stop
         import pr
         real(pr), intent(in out) :: X(:) !! Vector of variables \(X\)
         integer, intent(in out) :: ns !! Position of specified variable
         real(pr), intent(in out) :: S !! Specification variable value
         real(pr), intent(in out) :: dS !! Step of specification in the method
         real(pr), intent(in out) ::  dXdS(:) !! \(\frac{dX}{dS}\)
         integer, intent(in) :: iterations !! Iterations needed to converge point
      end function continuation_stopper
   end interface

   abstract interface
      subroutine continuation_solver(&
         fun, iters, X, ns, S, max_iters, F, dF, dFdS, tol &
         )
         !! Solver to solve a point during numerical contination.
         import pr, continuation_function
         procedure(continuation_function) :: fun !! Function to solve
         integer,  intent(out) :: iters !! Number of iterations needed
         real(pr), intent(in out) :: X(:)  !! Variables vector
         integer, intent(in) :: ns !! Specification number
         real(pr), intent(in) :: S !! Specification value
         integer, intent(in) :: max_iters !! Maximum iterations
         real(pr), intent(out) :: F(:) !! Function values at solved point
         real(pr), intent(out) :: df(:, :) !! Jacobian values
         real(pr), intent(out) :: dfds(:) !! dFdS
         real(pr), intent(in) :: tol !! Solver tolerance
      end subroutine continuation_solver
   end interface

contains

   function continuation(&
      f, X0, ns0, S0, dS0, max_points, solver_tol, &
      update_specification, postprocess, solver, stop &
      ) result(XS)
      !! Numerical continuation of a function.
      !!
      !! Uses Algower method of numerical continuation to trace a line that
      !! solves a system of the kind:
      !!
      !! \[ F(X,S) = 0 \]
      !!
      !! Where \(X\) is the variables vector and \(S)\ is the value of the
      !! specification.
      !! The method works with by providing a good set of initial points to
      !! solve the system of equations with an extrapolation using the previous
      !! solved point information.
      procedure(continuation_function) :: f !! Function to trace
      real(pr), intent(in) :: X0(:) !! Initial point
      integer, intent(in) :: ns0 !! Initial specification
      real(pr), intent(in) :: S0 !! Initial specification value
      real(pr), intent(in) :: dS0 !! Initial \(\deltaS\)
      integer, intent(in) :: max_points !! Maximum number of points to trace
      real(pr), intent(in) :: solver_tol !! Point solver tolerance
      procedure(process), optional :: update_specification
      !! Procedure to select the new specification and define the next step
      !! \(\DeltaS)\, defaults to:
      !!
      !! ```fortran
      !! ns = maxloc(abs(dXdS), dim=1)
      !! dS = dXdS(ns)*dS
      !! dXdS = dXdS/dXdS(ns)
      !! dS = sign(minval(abs([0.05_pr, dS])), dS)
      !! ```
      procedure(process), optional :: postprocess
      !! Any kind of postprocess that could be done after defining the
      !! next step
      procedure(continuation_solver), optional :: solver
      !! Solver procedures, uses Newton-Raphson by default
      procedure(continuation_stopper), optional :: stop
      !! Stopping procedure
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
            call solver(&
               f, newton_its, X, ns, S, max_iters, fval, dF, dFdS, solver_tol&
               )
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

         if (present(stop)) then
            if (stop(X, ns, S, dS, dXdS, newton_its)) exit
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
      procedure(continuation_function) :: fun !! Function to solve
      integer,  intent(out) :: iters !! Number of iterations needed
      real(pr), intent(in out) :: X(:)  !! Variables vector
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      integer, intent(in) :: max_iters !! Maximum iterations
      real(pr), intent(out) :: F(:) !! Function values at solved point
      real(pr), intent(out) :: df(:, :) !! Jacobian values
      real(pr), intent(out) :: dfds(:) !! dFdS
      real(pr), optional, intent(in) :: tol

      real(pr) :: dX(size(X)), solve_tol

      solve_tol = optval(tol, 1.e-5_pr)

      dX = 20
      F = 500
      newton: do iters = 1, max_iters
         ! Converged point
         if (maxval(abs(dx/x)) < solve_tol .or. maxval(abs(F)) < solve_tol) exit newton

         call fun(X, ns, S, F, dF, dFdS)
         if (maxval(abs(F)) < solve_tol) exit

         dX = solve_system(dF, -F)

         do while(maxval(abs(dx)) > 1)
            dX = dX/2
         end do

         X = X + dX
      end do newton
   end subroutine full_newton
end module yaeos__math_continuation
