module test_math_utilities_callbacks
   !! Callbacks with analytic solutions for the generic math utilities.
   use yaeos__constants, only: pr
   use yaeos__math_linalg, only: solve_system
   use yaeos__math_continuation, only: continuation_function
   implicit none

   integer :: solver_calls = 0
   integer :: postprocess_calls = 0
   integer :: stop_calls = 0
   integer :: stop_after = 0
   real(pr) :: fixed_dS = 0.0_pr

contains

   function vec_fun(x) result(f)
      !! f(x) = [x1^2 x2, sin(x1) + x2^3]
      real(pr), intent(in) :: x(:)
      real(pr) :: f(size(x))
      f(1) = x(1)**2 * x(2)
      f(2) = sin(x(1)) + x(2)**3
   end function vec_fun

   subroutine sqrt2_fun(x, f, df)
      !! f(x) = x^2 - 2, root sqrt(2)
      real(pr), intent(in) :: x
      real(pr), intent(out) :: f
      real(pr), intent(out) :: df
      f = x**2 - 2.0_pr
      df = 2.0_pr * x
   end subroutine sqrt2_fun

   subroutine hybrid_fun(n, x, fvec, iflag)
      !! Square system with root [1, 2]:
      !!   x1^2 + x2^2 - 5 = 0,  x1 x2 - 2 = 0
      integer, intent(in) :: n
      real(pr), intent(in) :: x(n)
      real(pr), intent(out) :: fvec(n)
      integer, intent(in out) :: iflag
      fvec(1) = x(1)**2 + x(2)**2 - 5.0_pr
      fvec(2) = x(1) * x(2) - 2.0_pr
   end subroutine hybrid_fun

   subroutine circle(X, ns, S, F, dF, dFdS)
      !! Unit circle with specification X(ns) = S
      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: dF(:, :)
      real(pr), intent(out) :: dFdS(:)
      F(1) = X(1)**2 + X(2)**2 - 1.0_pr
      F(2) = X(ns) - S
      dF = 0.0_pr
      dF(1, :) = [2.0_pr*X(1), 2.0_pr*X(2)]
      dF(2, ns) = 1.0_pr
      dFdS = 0.0_pr
      dFdS(2) = -1.0_pr
   end subroutine circle

   subroutine plain_newton(fun, iters, X, ns, S, dS, dXdS, point, max_iters, F, dF, dFdS, tol)
      !! Undamped Newton corrector used as a user-supplied `solver`.
      procedure(continuation_function) :: fun
      integer, intent(out) :: iters
      real(pr), intent(in out) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in) :: dS
      real(pr), intent(in) :: dXdS(:)
      integer, intent(in) :: point
      integer, intent(in) :: max_iters
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: dF(:, :)
      real(pr), intent(out) :: dFdS(:)
      real(pr), intent(in) :: tol
      solver_calls = solver_calls + 1
      do iters = 1, max_iters
         call fun(X, ns, S, F, dF, dFdS)
         if (maxval(abs(F)) < tol) exit
         X = X + solve_system(dF, -F)
      end do
   end subroutine plain_newton

   subroutine keep_spec(X, ns, S, dS, dXdS, iterations)
      !! Keep specifying X(2) with a constant step `fixed_dS`.
      real(pr), intent(in out) :: X(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: dS
      real(pr), intent(in out) :: dXdS(:)
      integer, intent(in) :: iterations
      ns = 2
      dS = fixed_dS
   end subroutine keep_spec

   subroutine count_points(X, ns, S, dS, dXdS, iterations)
      real(pr), intent(in out) :: X(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: dS
      real(pr), intent(in out) :: dXdS(:)
      integer, intent(in) :: iterations
      postprocess_calls = postprocess_calls + 1
   end subroutine count_points

   logical function stop_counter(X, ns, S, dS, dXdS, iterations)
      real(pr), intent(in out) :: X(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: dS
      real(pr), intent(in out) :: dXdS(:)
      integer, intent(in) :: iterations
      stop_calls = stop_calls + 1
      stop_counter = stop_calls >= stop_after
   end function stop_counter
end module test_math_utilities_callbacks

program test_math_utilities
   use yaeos__constants, only: pr
   use yaeos__math, only: central_numdiff, newton, sq_error, dx_to_dn, &
      interpol, intersect_one_line, Point, powel_hybrid
   use yaeos__math_continuation, only: continuation
   use auxiliar_functions, only: allclose, rel_error
   use test_math_utilities_callbacks
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Math utilities")

   call test_central_numdiff()
   call test_newton_1d()
   call test_small_helpers()
   call test_allclose()
   call test_intersections()
   call test_powell_hybrid()
   call test_continuation_callbacks()

   write(*, *) ""

contains

   subroutine test_central_numdiff()
      real(pr) :: x(2)
      real(pr) :: jac(2, 2)
      real(pr) :: exact(2, 2)

      x = [0.7_pr, 1.3_pr]
      jac = central_numdiff(vec_fun, x, 1e-5_pr)
      exact(1, :) = [2*x(1)*x(2), x(1)**2]
      exact(2, :) = [cos(x(1)), 3*x(2)**2]
      call assert(maxval(abs(jac - exact)) < 1e-8_pr, &
         "central_numdiff: matches analytic Jacobian to O(eps^2)")
   end subroutine test_central_numdiff

   subroutine test_newton_1d()
      real(pr) :: x
      logical :: failed

      x = 1.0_pr
      call newton(sqrt2_fun, x, 1e-12_pr, 50, failed)
      call assert(.not. failed, "newton_1d: converges for x^2 - 2")
      call assert(abs(x - sqrt(2.0_pr)) < 1e-10_pr, "newton_1d: root sqrt(2)")

      ! Large initial step is halved until |step| <= |x|/2
      x = 0.1_pr
      call newton(sqrt2_fun, x, 1e-12_pr, 100, failed)
      call assert(.not. failed, "newton_1d: converges from x0 = 0.1 with step halving")
      call assert(abs(x - sqrt(2.0_pr)) < 1e-10_pr, "newton_1d: root sqrt(2) from x0 = 0.1")

      x = 1.0_pr
      call newton(sqrt2_fun, x, 1e-12_pr, 1, failed)
      call assert(failed, "newton_1d: one iteration is not enough, failed is set")
      call assert(abs(x - 1.5_pr) < 1e-14_pr, "newton_1d: one Newton step from 1 gives 1.5")
   end subroutine test_newton_1d

   subroutine test_small_helpers()
      real(pr) :: x(3)
      real(pr) :: dx(3)
      real(pr) :: dn(3)

      call assert(abs(sq_error(2.0_pr, 1.0_pr) - 0.25_pr) < 1e-15_pr, "sq_error: ((2-1)/2)^2 = 0.25")
      call assert(sq_error(3.0_pr, 3.0_pr) == 0.0_pr, "sq_error: zero for exact prediction")

      ! dx_to_dn: dn = dx - sum(x dx); the result is orthogonal to x
      x = [0.2_pr, 0.3_pr, 0.5_pr]
      dx = [1.0_pr, -2.0_pr, 4.0_pr]
      dn = dx_to_dn(x, dx)
      call assert(maxval(abs(dn - (dx - 1.6_pr))) < 1e-14_pr, "dx_to_dn: dn = dx - sum(x*dx)")
      call assert(abs(sum(x*dn)) < 1e-14_pr, "dx_to_dn: sum(x*dn) = 0 for normalized x")

      call assert(abs(interpol(1.0_pr, 3.0_pr, 10.0_pr, 20.0_pr, 2.0_pr) - 15.0_pr) < 1e-14_pr, &
         "interpol: midpoint of a line")
      call assert(abs(interpol(1.0_pr, 3.0_pr, 10.0_pr, 20.0_pr, 5.0_pr) - 30.0_pr) < 1e-13_pr, &
         "interpol: linear extrapolation")
   end subroutine test_small_helpers

   subroutine test_allclose()
      call assert(abs(rel_error(2.0_pr, 2.5_pr) - 0.25_pr) < 1e-15_pr, "rel_error: |2-2.5|/2")
      call assert(allclose([1.0_pr, 2.0_pr], [1.0_pr, 2.0_pr + 1e-10_pr], 1e-8_pr), &
         "allclose: true within tolerance")
      call assert(.not. allclose([1.0_pr, 2.0_pr], [1.0_pr, 2.1_pr], 1e-3_pr), &
         "allclose: false when relative error 0.05 exceeds rtol")
   end subroutine test_allclose

   subroutine test_intersections()
      type(Point), allocatable :: inter(:)
      real(pr), allocatable :: lx(:)
      real(pr), allocatable :: ly(:)

      ! A horizontal segment crossed twice by the returning zigzag at x = 17 and 13
      lx = [0.0_pr, 20.0_pr, 18.0_pr, 16.0_pr, 14.0_pr, 12.0_pr]
      ly = [0.0_pr, 0.0_pr, 2.0_pr, -2.0_pr, -2.0_pr, 2.0_pr]
      inter = intersect_one_line(lx, ly)
      call assert(size(inter) == 2, "intersections: two crossings found")
      call assert(abs(inter(1)%x - 17.0_pr) < 1e-12_pr .and. abs(inter(1)%y) < 1e-12_pr, &
         "intersections: first crossing at (17, 0)")
      call assert(abs(inter(2)%x - 13.0_pr) < 1e-12_pr .and. abs(inter(2)%y) < 1e-12_pr, &
         "intersections: second crossing at (13, 0)")
      call assert(inter(1)%i == 2 .and. inter(1)%j == 3, "intersections: segment indices of first crossing")

      ! More than three crossings: the result is discarded (empty)
      lx = [0.0_pr, 20.0_pr, 18.0_pr, 16.0_pr, 14.0_pr, 12.0_pr, 10.0_pr, 8.0_pr]
      ly = [0.0_pr, 0.0_pr, 2.0_pr, -2.0_pr, 2.0_pr, -2.0_pr, 2.0_pr, -2.0_pr]
      inter = intersect_one_line(lx, ly)
      call assert(size(inter) == 0, "intersections: more than three crossings returns none")
   end subroutine test_intersections

   subroutine test_powell_hybrid()
      real(pr) :: x(2)
      real(pr) :: f(2)
      integer :: info

      x = [0.5_pr, 2.5_pr]
      call powel_hybrid(hybrid_fun, 1e-12_pr, x, f, info)
      call assert(info == 1, "powel_hybrid: MINPACK info = 1 (relative error below tol)")
      call assert(maxval(abs(x - [1.0_pr, 2.0_pr])) < 1e-9_pr, "powel_hybrid: root [1, 2]")
      call assert(maxval(abs(f)) < 1e-9_pr, "powel_hybrid: residual at the root is ~0")
   end subroutine test_powell_hybrid

   subroutine test_continuation_callbacks()
      integer, parameter :: max_points = 20
      real(pr) :: XS(max_points, 2)
      integer :: i

      ! Custom corrector + custom spec update + postprocess + stop
      solver_calls = 0
      postprocess_calls = 0
      stop_calls = 0
      stop_after = 6
      fixed_dS = -0.05_pr
      XS = continuation(circle, X0=[0.6_pr, 0.8_pr], ns0=2, S0=0.8_pr, dS0=-0.05_pr, &
         max_points=max_points, solver_tol=1e-12_pr, &
         update_specification=keep_spec, postprocess=count_points, &
         solver=plain_newton, stop=stop_counter)

      call assert(solver_calls == 6, "continuation: user solver used for every point")
      call assert(postprocess_calls == 6, "continuation: postprocess called once per point")
      call assert(stop_calls == 6, "continuation: stop evaluated once per point until true")
      do i = 1, 6
         call assert(abs(XS(i, 2) - (0.8_pr - 0.05_pr*(i - 1))) < 1e-10_pr, &
            "continuation: constant step in X(2) from update_specification")
         call assert(abs(XS(i, 1)**2 + XS(i, 2)**2 - 1.0_pr) < 1e-10_pr, &
            "continuation: every stored point lies on the unit circle")
         call assert(XS(i, 1) > 0.0_pr, "continuation: stays on the x1 > 0 branch")
      end do
      call assert(all(XS(7:, :) == 0.0_pr), "continuation: rows after stop remain zero")

      ! Default solver and specification, stopped by the user after 4 points
      stop_calls = 0
      stop_after = 4
      XS = continuation(circle, X0=[0.6_pr, 0.8_pr], ns0=2, S0=0.8_pr, dS0=-0.05_pr, &
         max_points=max_points, solver_tol=1e-10_pr, stop=stop_counter)
      do i = 1, 4
         call assert(abs(XS(i, 1)**2 + XS(i, 2)**2 - 1.0_pr) < 1e-8_pr, &
            "continuation (default solver): stored points lie on the unit circle")
      end do
      call assert(all(XS(5:, :) == 0.0_pr), "continuation (default solver): rows after stop remain zero")
   end subroutine test_continuation_callbacks
end program test_math_utilities
