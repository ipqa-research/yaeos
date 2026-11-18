module yaeos__math_continuation
   !! # Numerical continuation
   !! Predictor-corrector continuation of nonlinear solution curves.
   !!
   !! # Description
   !! The module defines callback interfaces and a continuation driver for
   !! tracing solutions of \(F(X,S)=0\) as a scalar specification \(S\)
   !! changes. At each accepted point the routine solves for the local tangent
   !! \(dX/dS\), predicts the next state, and corrects it with a nonlinear
   !! solver. The caller may customize the specification update, post-processing,
   !! point solver, and stopping condition.
   !!
   !! # Examples
   !! ```fortran
   !! curve = continuation(fun, X0, ns0, S0, dS0, &
   !!    max_points=100, solver_tol=1.0e-8_pr)
   !! ```
   use yaeos__constants, only: pr
   use yaeos__math_linalg, only: solve_system

   implicit none

   type :: ContinuationVariable
      !! # ContinuationVariable
      !! State of the unknown vector and specification at one curve point.
      !!
      !! # Description
      !! Records which component is specified, that component's value and the
      !! proposed step for the next point.
      !!
      !! # Examples
      !! ```fortran
      !! point%X = [P, T, V]
      !! point%ns = 1
      !! point%S = P
      !! point%dS = 0.1_pr
      !! ```
      real(pr), allocatable :: X(:) !! Unknown vector at the point.
      integer :: ns !! Index of the component of X used as the specification.
      real(pr) :: S !! Value of the specified variable.
      real(pr) :: dS !! Signed step proposed for the next continuation point.
   end type ContinuationVariable

   abstract interface
      subroutine continuation_function(X, ns, S, F, dF, dFdS)
         !! # continuation_function
         !! Evaluate the nonlinear residual and its local derivatives.
         !!
         !! # Description
         !! The equations are \(F(X,S)=0\). Return the Jacobian with respect
         !! to \(X\) and the derivative with respect to the scalar \(S\).
         import pr
         real(pr), intent(in) :: X(:) !! Current unknown vector.
         integer, intent(in) :: ns !! Index of the specified entry in X.
         real(pr), intent(in) :: S !! Current specification value.
         real(pr), intent(out) :: F(:) !! Residual vector, zero at a solution.
         real(pr), intent(out) :: dF(:, :) !! Jacobian \(\partial F/\partial X\).
         real(pr), intent(out) :: dFdS(:) !! Derivative \(\partial F/\partial S\).
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
         !! # continuation_stopper
         !!
         !! # Description
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
         fun, iters, X, ns, S, dS, dXdS, point, max_iters, F, dF, dFdS, tol &
         )
         !! # continuation_solver
         !!
         !! # Description
         !! Solver to solve a point during numerical continuation.
         import pr, continuation_function
         procedure(continuation_function) :: fun !! Function to solve
         integer,  intent(out) :: iters !! Number of iterations needed
         real(pr), intent(in out) :: X(:) !! Variables vector
         integer, intent(in) :: ns !! Specification number
         real(pr), intent(in) :: S !! Specification value
         real(pr), intent(in) :: dS !! Delta spec
         real(pr), intent(in) :: dXdS(:) !! Tangent \(dX/dS\) used for the predictor
         integer, intent(in) :: point !! Point number
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
      !! # continuation
      !! Trace a one-dimensional solution manifold of a nonlinear system.
      !!
      !! # Description
      !! At each point the tangent is obtained from
      !!
      !! \[ \frac{dX}{dS} = -\left(\frac{\partial F}{\partial X}\right)^{-1}
      !!       \frac{\partial F}{\partial S}. \]
      !!
      !! The tangent predicts the next point; `solver` (or the default Newton
      !! corrector) then solves the residual equations there. The default
      !! specification update selects the component of the tangent with the
      !! largest magnitude, rescales the tangent to it and converts the step to
      !! that variable (no step-size limit is applied). Tracing stops when the
      !! corrector reaches its iteration limit (500), the residual is NaN,
      !! `dS` becomes zero, or `stop` returns true. The returned array has `max_points` rows; rows not reached
      !! before termination remain zero.
      !!
      !! # Examples
      !! ```fortran
      !! curve = continuation(fun, X0, ns0, S0, dS0, &
      !!    max_points=100, solver_tol=1.0e-8_pr)
      !! ```
      procedure(continuation_function) :: f !! Residual and derivative callback for \(F(X,S)=0\).
      real(pr), intent(in) :: X0(:) !! Initial converged unknown vector.
      integer, intent(in) :: ns0 !! Index of the initially specified variable in `X0`.
      real(pr), intent(in) :: S0 !! Initial value of the specified variable.
      real(pr), intent(in) :: dS0 !! Signed first step in the specification.
      integer, intent(in) :: max_points !! Maximum number of continuation points to store.
      real(pr), intent(in) :: solver_tol !! Residual tolerance for the point corrector.
      procedure(process), optional :: update_specification
      !! Optional callback to select the next specification and step. If absent,
      !! the tangent component of largest magnitude is selected and the step is
      !! scaled using that component:
      !!
      !! ```fortran
      !! ns = maxloc(abs(dXdS), dim=1)
      !! dS = dXdS(ns)*dS
      !! dXdS = dXdS/dXdS(ns)
      !! ```
      procedure(process), optional :: postprocess
      !! Optional callback to modify the state or step after specification
      !! selection, for example to enforce bounds or change step direction.
      procedure(continuation_solver), optional :: solver
      !! Optional point corrector. The built-in Newton solver is used by default.
      procedure(continuation_stopper), optional :: stop
      !! Optional callback evaluated after each accepted point; true terminates
      !! continuation.
      real(pr) :: XS(max_points, size(X0)) !! Stored solution vectors, one per row.

      real(pr) :: X(size(X0))
      real(pr) :: S
      real(pr) :: fval(size(X0))
      real(pr) :: dF(size(X0), size(X0))
      real(pr) :: dFdS(size(X0))
      real(pr) :: dXdS(size(X0))
      integer :: ns
      real(pr) :: dS

      integer :: i !! Loop/index variable.
      integer :: newton_its

      integer :: max_iters = 500

      X = X0
      ns = ns0
      dS = dS0
      S = S0

      XS = 0

      do i=1,max_points

         if (present(solver)) then
            call solver(&
               f, newton_its, X, ns, S, dS, dXdS, i, max_iters, &
               fval, dF, dFdS, solver_tol&
               )
         else
            call full_newton(&
               f, newton_its, X, ns, S, dS, dXdS, i, max_iters, &
               fval, dF, dFdS, solver_tol &
               )
         end if
         if (newton_its >= max_iters .or. dS == 0._pr .or. any(isnan(Fval)))  exit

         XS(i, :) = X

         dXdS = solve_system(dF, -dFdS)

         if (present(update_specification)) then
            call update_specification(X, ns, S, dS, dXdS, newton_its)
         else
            ns = maxloc(abs(dXdS), dim=1)
            dS = dXdS(ns)*dS
            dXdS = dXdS/dXdS(ns)
         end if

         if (present(postprocess)) then
            call postprocess(X, ns, S, dS, dXdS, newton_its)
         end if

         if (present(stop)) then
            if (stop(X, ns, S, dS, dXdS, newton_its)) exit
         end if

         if (dS == 0) exit

         X = X + dXdS * dS
         S = X(ns)
      end do
   end function continuation

   subroutine full_newton(&
      fun, iters, X, ns, S, dS, dXdS, point, max_iters, F, dF, dFdS, tol &
      )
      !! # full_newton
      !!
      !! # Description
      !! Correct one predicted continuation point with Newton's method.
      !! The correction solves \(F(X,S)=0\) at fixed specification value,
      !! halving the Newton step until \(\max|\Delta X_i| \le 0.1\). It stops
      !! when \(\max|\Delta X_i/X_i| < 10^{-5}\) or \(\max|F_i| <\) `tol`.
      !! If iteration 10 is reached (for points after the first), \(X\) is reset
      !! to the previous point plus 10% of the predictor step. On failure
      !! `iters` ends at `max_iters + 1`.
      use yaeos__auxiliar, only: optval
      use yaeos__math_linalg, only: solve_system
      procedure(continuation_function) :: fun !! Residual/Jacobian callback.
      integer,  intent(out) :: iters !! Newton iterations used for this point.
      real(pr), intent(in out) :: X(:) !! Predicted vector; corrected solution on return.
      integer, intent(in) :: ns !! Index of the specified variable.
      real(pr), intent(in) :: S !! Fixed specification value for this correction.
      real(pr), intent(in) :: dS !! Step used to form the predictor.
      real(pr), intent(in) :: dXdS(:) !! Tangent used by the predictor and recovery logic.
      integer, intent(in) :: point !! One-based continuation point number.
      integer, intent(in) :: max_iters !! Maximum Newton iterations for correction.
      real(pr), intent(out) :: F(:) !! Residual at the corrected point.
      real(pr), intent(out) :: df(:, :) !! Jacobian \(\partial F/\partial X\).
      real(pr), intent(out) :: dfds(:) !! Derivative \(\partial F/\partial S\).
      real(pr), intent(in) :: tol !! Residual (max-norm) convergence tolerance.

      real(pr) :: X0(size(X))

      real(pr) :: dX(size(X))
      real(pr) :: solve_tol

      solve_tol = tol

      dX = 20
      F = 500
      X0 = X
      newton: do iters = 1, max_iters
         call fun(X, ns, S, F, dF, dFdS)
         dX = solve_system(dF, -F)

         ! Converged point
         if (maxval(abs(dx/x)) < 1e-5_pr .or. maxval(abs(F)) < tol) exit newton

         ! Fix the step
         do while(maxval(abs(dx)) > 0.1)
            dX = dX/2
         end do

         if (iters == 10 .and. point > 1) then
            dX = X0 - 0.9*dXdS*dS - X
         end if

         X = X + dX
      end do newton
   end subroutine full_newton

   ! subroutine levenberg_marquardt(&
   !    fun, iters, X, ns, S, dS, dXdS, point, max_iters, F, dF, dFdS, tol &
   !    )
   !    use minpack_module, only: lmdif1
   !    use yaeos__auxiliar, only: optval
   !    use yaeos__math_linalg, only: solve_system
   !    procedure(continuation_function) :: fun !! Function to solve
   !    integer,  intent(out) :: iters !! Number of iterations needed
   !    real(pr), intent(in out) :: X(:)  !! Variables vector
   !    integer, intent(in) :: ns
   !    real(pr), intent(in) :: S
   !    real(pr), intent(in) :: dS
   !    real(pr), intent(in) :: dXdS(:)
   !    integer, intent(in) :: point
   !    integer, intent(in) :: max_iters !! Maximum iterations
   !    real(pr), intent(out) :: F(:) !! Function values at solved point
   !    real(pr), intent(out) :: df(:, :) !! Jacobian values
   !    real(pr), intent(out) :: dfds(:) !! dFdS
   !    real(pr), intent(in) :: tol

   !    integer :: m, n, info, iwa(size(x))
   !    integer :: lwa
   !    real(pr) :: wa(size(F) * size(x) + 5*size(x) + size(f))

   !    m = size(F)
   !    n = size(x)
   !    lwa = size(F) * size(x)+5*size(x)+size(f)
   !    call lmdif1(fcn, m, n, x, F, tol, Info, Iwa, Wa, Lwa)
   !    contains
   !    subroutine fcn(m, n, xx, fvec, iflag)
   !       integer, intent(in) :: m, n
   !       real(pr), intent(in) :: xx(n)
   !       real(pr), intent(out) :: fvec(m)
   !       integer, intent(in out) :: iflag
   !       call fun(xx, ns, S, fvec, dF, dFdS)
   !    end subroutine
   ! end subroutine
end module yaeos__math_continuation
