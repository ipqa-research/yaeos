module test_newton_solver_systems
   !! Residual/Jacobian callbacks with known roots used to exercise the
   !! safeguarded Newton-LM solver in `yaeos__newton_solver`.
   use yaeos, only: pr
   use iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none

   integer :: n_calls = 0

contains

   subroutine linear_2x2(X, F, J)
      !! F = A x - b with A = [[2, 1], [1, 3]], b = [3, 5]  ->  root x = [0.8, 1.4]
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      n_calls = n_calls + 1
      J(1, :) = [2.0_pr, 1.0_pr]
      J(2, :) = [1.0_pr, 3.0_pr]
      F = matmul(J, X) - [3.0_pr, 5.0_pr]
   end subroutine linear_2x2

   subroutine cubic_1d(X, F, J)
      !! F = x^3: Newton step is exactly dX = -x/3, so x_k = x_0 (2/3)^k
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      F(1) = X(1)**3
      J(1, 1) = 3.0_pr * X(1)**2
   end subroutine cubic_1d

   subroutine atan_1d(X, F, J)
      !! F = atan(x): pure Newton overshoots from |x0| > 1.39, root x = 0
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      F(1) = atan(X(1))
      J(1, 1) = 1.0_pr / (1.0_pr + X(1)**2)
   end subroutine atan_1d

   subroutine log_1d(X, F, J)
      !! F = log(x) (NaN for x <= 0), root x = 1
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      if (X(1) > 0.0_pr) then
         F(1) = log(X(1))
      else
         F(1) = ieee_value(1.0_pr, ieee_quiet_nan)
      end if
      J(1, 1) = 1.0_pr / X(1)
   end subroutine log_1d

   subroutine nan_start(X, F, J)
      !! Residual is NaN everywhere
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      F = ieee_value(1.0_pr, ieee_quiet_nan)
      J = 0.0_pr
      J(1, 1) = 1.0_pr
      J(2, 2) = 1.0_pr
      F(2) = X(2)
   end subroutine nan_start

   subroutine rank_deficient(X, F, J)
      !! F = [x1 + x2 - 2, 2 (x1 + x2 - 2)]: singular Jacobian everywhere,
      !! the solution set is the line x1 + x2 = 2.
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      F(1) = X(1) + X(2) - 2.0_pr
      F(2) = 2.0_pr * F(1)
      J(1, :) = [1.0_pr, 1.0_pr]
      J(2, :) = [2.0_pr, 2.0_pr]
   end subroutine rank_deficient

   subroutine huge_rank_deficient(X, F, J)
      !! F = s [g, g] with g = x1 + x2 - 2 and s = 2^40. J = s [[1, 1], [1, 1]]
      !! (powers of two, so LU arithmetic is exact): a tiny LM shift is
      !! absorbed by J^T J and the regularized matrix stays exactly singular.
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      real(pr), parameter :: scale = 2.0_pr**40
      F = scale * (X(1) + X(2) - 2.0_pr)
      J = scale
   end subroutine huge_rank_deficient

   subroutine rosenbrock_root(X, F, J)
      !! F = [10 (x2 - x1^2), 1 - x1], root x = [1, 1]
      real(pr), intent(in) :: X(:)
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: J(:, :)
      F(1) = 10.0_pr * (X(2) - X(1)**2)
      F(2) = 1.0_pr - X(1)
      J(1, :) = [-20.0_pr * X(1), 10.0_pr]
      J(2, :) = [-1.0_pr, 0.0_pr]
   end subroutine rosenbrock_root
end module test_newton_solver_systems

program test_newton_solver
   use yaeos, only: pr
   use yaeos__newton_solver, only: newton, newton_settings, newton_result, &
      newton_status_msg, NEWTON_SUCCESS, NEWTON_MAX_ITS, NEWTON_NAN_INF, &
      NEWTON_SINGULAR, NEWTON_STAGNATION, NEWTON_LINE_SEARCH_FAIL
   use test_newton_solver_systems
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Safeguarded Newton-LM solver")

   call test_linear_one_step()
   call test_cubic_max_iterations()
   call test_cubic_stagnation()
   call test_nan_initial_point()
   call test_backtracking_atan()
   call test_ls_failure_promotes_lm()
   call test_t_min_failure()
   call test_line_search_ceiling()
   call test_nan_during_backtracking()
   call test_singular_jacobian_lm()
   call test_ill_conditioned_promotion()
   call test_forced_lm_from_start()
   call test_singular_regularized_matrix()
   call test_status_messages()

   write(*, *) ""

contains

   subroutine test_linear_one_step()
      !! A linear system is solved exactly by one Newton step; history and
      !! verbose printing are exercised.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(2)

      x = [10.0_pr, -7.0_pr]
      s%save_history = .true.
      s%verbosity = 2
      n_calls = 0
      call newton(linear_2x2, x, s, r)

      call assert(r%status == NEWTON_SUCCESS, "linear: status is NEWTON_SUCCESS")
      call assert(maxval(abs(x - [0.8_pr, 1.4_pr])) < 1e-12_pr, "linear: exact root [0.8, 1.4]")
      call assert(r%f_norm < s%atol, "linear: final residual norm below atol")
      call assert(r%lambda == 0.0_pr, "linear: pure Newton step, lambda stays zero")
      call assert(r%cond_est >= 1.0_pr .and. r%cond_est < 10.0_pr, &
         "linear: condition estimate of a well conditioned 2x2 matrix")
      call assert(allocated(r%f_history), "linear: history allocated with save_history")
      call assert(size(r%f_history) == r%iterations + 1, "linear: history has iterations+1 entries")
      ! Only the first two entries (x0 and after the single step) are checked:
      ! the trailing entry on a convergence exit is never written by the solver.
      call assert(abs(r%f_history(1) - maxval(abs([2*10.0_pr - 7.0_pr - 3, 10.0_pr - 21.0_pr - 5]))) < 1e-12_pr, &
         "linear: history starts with ||F(x0)||_inf")
      call assert(r%f_history(2) < 1e-12_pr, "linear: second history entry is ~0 after one exact step")
      call assert(n_calls == 2, "linear: one initial evaluation plus one trial evaluation")

      ! Default settings and no result: same root
      x = [10.0_pr, -7.0_pr]
      call newton(linear_2x2, x)
      call assert(maxval(abs(x - [0.8_pr, 1.4_pr])) < 1e-12_pr, "linear: same root with all defaults")

      ! Starting at the root converges with zero iterations of work
      x = [0.8_pr, 1.4_pr]
      call newton(linear_2x2, x, result=r)
      call assert(r%status == NEWTON_SUCCESS, "linear: starting at the root converges")
      call assert(all(x == [0.8_pr, 1.4_pr]), "linear: starting at the root, x is not moved")
      call assert(.not. allocated(r%f_history), "linear: no history without save_history")
   end subroutine test_linear_one_step

   subroutine test_cubic_max_iterations()
      !! On F = x^3 Newton converges linearly, x_k = x0 (2/3)^k
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)
      integer :: k

      x = [1.0_pr]
      s%max_its = 10
      s%save_history = .true.
      s%verbosity = 1
      call newton(cubic_1d, x, s, r)

      call assert(r%status == NEWTON_MAX_ITS, "cubic: status is NEWTON_MAX_ITS after max_its")
      call assert(r%iterations == 10, "cubic: iterations equals max_its")
      call assert(abs(x(1) - (2.0_pr/3.0_pr)**10) < 1e-12_pr, "cubic: x_10 = (2/3)^10")
      call assert(abs(r%dx_norm - (2.0_pr/3.0_pr)**9 / 3) < 1e-12_pr, "cubic: last step is x_9/3")
      call assert(size(r%f_history) == 11, "cubic: history of max_its + 1 entries")
      do k = 0, 10
         call assert(abs(r%f_history(k + 1) - (2.0_pr/3.0_pr)**(3*k)) < 1e-12_pr, &
            "cubic: history entry k equals |x_k|^3")
      end do
      call assert(abs(r%f_norm - (2.0_pr/3.0_pr)**30) < 1e-12_pr, "cubic: final f_norm is |x_10|^3")
   end subroutine test_cubic_max_iterations

   subroutine test_cubic_stagnation()
      !! With a generous stagnation tolerance every step after the first counts
      !! as stagnant, so the solver must stop at iteration stagnation_nits + 1.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)

      x = [1.0_pr]
      s%atol = 1e-300_pr
      s%rtol = 1e-300_pr
      s%stagnation_tol = 1.0_pr
      s%stagnation_nits = 3
      s%verbosity = 1
      call newton(cubic_1d, x, s, r)

      call assert(r%status == NEWTON_STAGNATION, "stagnation: status is NEWTON_STAGNATION")
      call assert(r%iterations == 4, "stagnation: detected at iteration stagnation_nits + 1")
      call assert(abs(x(1) - (2.0_pr/3.0_pr)**3) < 1e-12_pr, "stagnation: three accepted steps, x = (2/3)^3")
      call assert(newton_status_msg(r%status) == "Stagnation: no progress detected", &
         "stagnation: status message")
   end subroutine test_cubic_stagnation

   subroutine test_nan_initial_point()
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(2)

      x = [1.0_pr, 2.0_pr]
      s%save_history = .true.
      call newton(nan_start, x, s, r)
      call assert(r%status == NEWTON_NAN_INF, "nan start: status is NEWTON_NAN_INF")
      call assert(r%iterations == 0, "nan start: no iterations performed")
      call assert(r%f_norm == huge(1.0_pr), "nan start: f_norm reported as huge")
      call assert(all(x == [1.0_pr, 2.0_pr]), "nan start: x left untouched")
      call assert(size(r%f_history) == 1, "nan start: history has a single entry")

      ! Without the optional result it must simply return
      x = [1.0_pr, 2.0_pr]
      call newton(nan_start, x)
      call assert(all(x == [1.0_pr, 2.0_pr]), "nan start without result: x left untouched")
   end subroutine test_nan_initial_point

   subroutine test_backtracking_atan()
      !! From x0 = 3 the full Newton step on atan(x) overshoots to |x| > 3;
      !! Armijo backtracking must still bring it to the root.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)
      real(pr) :: f(1)
      real(pr) :: jac(1, 1)

      x = [3.0_pr]
      call atan_1d(x, f, jac)
      call assert(abs(x(1) - f(1)/jac(1, 1)) > 3.0_pr, "atan: full Newton step from x0 = 3 overshoots")

      s%verbosity = 2
      call newton(atan_1d, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "atan: backtracking reaches convergence")
      call assert(abs(x(1)) < 1e-9_pr, "atan: root x = 0")
      call assert(r%lambda == 0.0_pr, "atan: Armijo succeeded so no LM was needed")
   end subroutine test_backtracking_atan

   subroutine test_ls_failure_promotes_lm()
      !! With a single backtracking step allowed the Armijo search fails, LM
      !! is activated (lambda > 0) and the solver still reaches the root.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)

      x = [3.0_pr]
      s%armijo_max_its = 1
      s%max_its = 200
      s%verbosity = 2
      call newton(atan_1d, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "ls failure: LM fallback converges")
      call assert(abs(x(1)) < 1e-9_pr, "ls failure: root x = 0 found with LM")
      call assert(r%lambda > 0.0_pr, "ls failure: lambda was activated")
   end subroutine test_ls_failure_promotes_lm

   subroutine test_t_min_failure()
      !! A minimum step length above armijo_tau makes the line search fail
      !! as soon as t would drop below t_min; LM then rescues the iteration.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)

      x = [3.0_pr]
      s%t_min = 0.6_pr
      s%max_its = 200
      call newton(atan_1d, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "t_min: LM fallback converges")
      call assert(abs(x(1)) < 1e-9_pr, "t_min: root x = 0")
      call assert(r%lambda > 0.0_pr, "t_min: lambda was activated by the line-search failure")
   end subroutine test_t_min_failure

   subroutine test_line_search_ceiling()
      !! lambda_max = lambda_min * lambda_up: the first line-search failure
      !! reaches the ceiling and the solver must give up at x0.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)
      real(pr) :: jscale

      x = [3.0_pr]
      s%armijo_max_its = 1
      s%lambda_max = 1.0e-5_pr
      s%verbosity = 1
      call newton(atan_1d, x, s, r)
      jscale = (1.0_pr / 10.0_pr)**2
      call assert(r%status == NEWTON_LINE_SEARCH_FAIL, "ceiling: status is NEWTON_LINE_SEARCH_FAIL")
      call assert(x(1) == 3.0_pr, "ceiling: last valid iterate (x0) returned")
      call assert(r%iterations == 1, "ceiling: failure on the first iteration")
      call assert(abs(r%lambda - jscale * s%lambda_max) < 1e-12_pr * jscale, &
         "ceiling: lambda equals jscale * lambda_max")
      call assert(abs(r%f_norm - atan(3.0_pr)) < 1e-14_pr, "ceiling: f_norm is |F(x0)|")
   end subroutine test_line_search_ceiling

   subroutine test_nan_during_backtracking()
      !! From x0 = 10 Newton on log(x) jumps to x < 0 (NaN); the line search
      !! must keep halving through NaN trial points.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(1)

      x = [10.0_pr]
      call newton(log_1d, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "log: converges through NaN trial points")
      call assert(abs(x(1) - 1.0_pr) < 1e-9_pr, "log: root x = 1")

      ! Single backtracking step that lands on NaN: line search fails, LM rescues
      x = [10.0_pr]
      s%armijo_max_its = 1
      s%max_its = 500
      call newton(log_1d, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "log with one backtrack: LM rescue converges")
      call assert(abs(x(1) - 1.0_pr) < 1e-9_pr, "log with one backtrack: root x = 1")
      call assert(r%lambda > 0.0_pr, "log with one backtrack: LM was activated")
   end subroutine test_nan_during_backtracking

   subroutine test_singular_jacobian_lm()
      !! Exactly singular J: dgetrf fails, the solver promotes to LM and finds
      !! a point on the solution line x1 + x2 = 2.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(2)

      x = [3.0_pr, 4.0_pr]
      s%max_its = 200
      call newton(rank_deficient, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "singular J: LM converges")
      call assert(abs(x(1) + x(2) - 2.0_pr) < 1e-8_pr, "singular J: point on x1 + x2 = 2")
      call assert(abs(x(1) - x(2) + 1.0_pr) < 1e-6_pr, &
         "singular J: LM steps are along J^T, so x1 - x2 is preserved (= -1)")
      call assert(r%lambda > 0.0_pr, "singular J: LM lambda active")
   end subroutine test_singular_jacobian_lm

   subroutine test_ill_conditioned_promotion()
      !! cond_max below any real condition number forces the LM branch on a
      !! nonsingular problem; the root must still be found.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(2)

      x = [-1.2_pr, 1.0_pr]
      s%cond_max = 0.5_pr
      s%max_its = 500
      call newton(rosenbrock_root, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "cond_max: LM-only iteration converges")
      call assert(maxval(abs(x - 1.0_pr)) < 1e-8_pr, "cond_max: root [1, 1]")
      call assert(r%lambda > 0.0_pr, "cond_max: lambda promoted by conditioning test")

      ! Reference: pure Newton with default settings gives the same root
      x = [-1.2_pr, 1.0_pr]
      call newton(rosenbrock_root, x, result=r)
      call assert(r%status == NEWTON_SUCCESS, "rosenbrock: default solver converges")
      call assert(maxval(abs(x - 1.0_pr)) < 1e-8_pr, "rosenbrock: default solver root [1, 1]")
   end subroutine test_ill_conditioned_promotion

   subroutine test_forced_lm_from_start()
      !! lambda0 > 0 means LM from the first iteration; lambda decays towards
      !! lambda0 but never below it.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(2)

      x = [10.0_pr, -7.0_pr]
      s%lambda0 = 1.0e-3_pr
      s%max_its = 200
      call newton(linear_2x2, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "lambda0: forced LM converges")
      call assert(maxval(abs(x - [0.8_pr, 1.4_pr])) < 1e-9_pr, "lambda0: root [0.8, 1.4]")
      call assert(r%lambda >= s%lambda0, "lambda0: lambda never decays below lambda0")
      call assert(r%iterations > 2, "lambda0: damped steps need more than one iteration")
   end subroutine test_forced_lm_from_start

   subroutine test_singular_regularized_matrix()
      !! With lambda_min = 1e-30 the LM shift is lost against J^T J ~ 2^81, so
      !! the first regularized matrix is exactly singular and the solver gets a
      !! zero step. A zero step is not convergence: the solver must ramp lambda
      !! and retry, reaching a real root on x1 + x2 = 2, or report
      !! NEWTON_SINGULAR when lambda cannot grow any further.
      type(newton_settings) :: s
      type(newton_result) :: r
      real(pr) :: x(2)

      ! Room to grow lambda: the ramped regularization becomes nonsingular.
      x = [3.0_pr, 4.0_pr]
      s%lambda_min = 1.0e-30_pr
      s%max_its = 500
      call newton(huge_rank_deficient, x, s, r)
      call assert(r%status == NEWTON_SUCCESS, "singular LM: recovers by ramping lambda")
      call assert(abs(x(1) + x(2) - 2.0_pr) < 1e-8_pr, "singular LM: success means a real root (x1 + x2 = 2)")

      ! No room to grow lambda: the matrix stays singular, so report it.
      x = [3.0_pr, 4.0_pr]
      s%lambda_min = 1.0e-30_pr
      s%lambda_max = 1.0e-30_pr
      call newton(huge_rank_deficient, x, s, r)
      call assert(r%status == NEWTON_SINGULAR, "singular LM at the lambda ceiling: NEWTON_SINGULAR")
      call assert(all(x == [3.0_pr, 4.0_pr]), "singular LM at the lambda ceiling: x unchanged")
      call assert(r%f_norm == 5.0_pr * 2.0_pr**40, "singular LM at the lambda ceiling: residual still |F(x0)|")
   end subroutine test_singular_regularized_matrix

   subroutine test_status_messages()
      call assert(newton_status_msg(NEWTON_SUCCESS) == "Converged", "msg: success")
      call assert(newton_status_msg(NEWTON_MAX_ITS) == "Maximum iterations reached", "msg: max its")
      call assert(newton_status_msg(NEWTON_NAN_INF) == "NaN / Inf detected in residual", "msg: nan/inf")
      call assert(newton_status_msg(NEWTON_SINGULAR) == "Singular Jacobian", "msg: singular")
      call assert(newton_status_msg(NEWTON_LINE_SEARCH_FAIL) == "Line search failed", "msg: line search")
      call assert(newton_status_msg(-42) == "Unknown status code", "msg: unknown code")
   end subroutine test_status_messages
end program test_newton_solver
