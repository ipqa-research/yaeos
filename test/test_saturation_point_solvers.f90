module test_saturation_point_solvers_callbacks
   !! Residual callbacks for `newton_solve`.
   use yaeos, only: pr
   implicit none
contains
   subroutine circle_line(x, f, J)
      !! \(x_1^2 + x_2^2 = 4\), \(x_1 = x_2\): root at \((\sqrt 2, \sqrt 2)\).
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: f(:)
      real(pr), intent(out) :: J(:, :)
      f(1) = x(1)**2 + x(2)**2 - 4
      f(2) = x(1) - x(2)
      J(1, :) = [2*x(1), 2*x(2)]
      J(2, :) = [1.0_pr, -1.0_pr]
   end subroutine circle_line

   subroutine arctangent(x, f, J)
      !! \(\arctan x = 0\). Undamped Newton diverges from \(|x_0| > 1.39\).
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: f(:)
      real(pr), intent(out) :: J(:, :)
      f(1) = atan(x(1))
      J(1, 1) = 1/(1 + x(1)**2)
   end subroutine arctangent
end module test_saturation_point_solvers_callbacks

program test_saturation_point_solvers
   !! Low-level saturation solvers of `src/models/solvers/saturation_point.f90`:
   !! the damped Newton solver, the log-volume saturation system
   !! (`saturation_F`, residuals and analytical Jacobian), its Newton solver
   !! (`solve_VxVyT`), and the stable-root branch of `saturation_TP`.
   !! The reference bubble point comes from `saturation_pressure`.
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, saturation_pressure
   use yaeos__m_s_sp, only: saturation_F, saturation_TP, solve_VxVyT
   use nonlinear_solvers, only: newton_solve
   use test_saturation_point_solvers_callbacks, only: circle_line, arctangent
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: bub
   real(pr) :: z(nc)
   real(pr) :: T

   write(*, *) test_title("Saturation point solvers")

   call test_newton_solve()

   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]
   T = 300.0_pr
   bub = saturation_pressure(model, z, T, kind="bubble")
   call assert(bub%kind == "bubble", "reference bubble point converges")

   call test_saturation_F()
   call test_solve_VxVyT()
   call test_saturation_TP_stable_roots()

contains

   subroutine test_newton_solve()
      real(pr) :: x2(2)
      real(pr) :: x1(1)
      integer :: its
      integer :: info

      x2 = [3.0_pr, 1.0_pr]
      call newton_solve(circle_line, x2, tol=1e-12_pr, max_iter=50, its=its, info=info)
      call assert(info == 0, "newton_solve: circle-line system converges")
      call assert(maxval(abs(x2 - sqrt(2.0_pr))) < 1e-10_pr, "newton_solve: root at (sqrt 2, sqrt 2)")

      ! Backtracking is needed: plain Newton on atan(x) diverges from x = 3
      x1 = [3.0_pr]
      call newton_solve(arctangent, x1, tol=1e-12_pr, max_iter=100, its=its, info=info)
      call assert(info == 0, "newton_solve: atan(x) = 0 converges from x0 = 3 with backtracking")
      call assert(abs(x1(1)) < 1e-12_pr, "newton_solve: atan root at zero")

      ! Iteration limit reached
      x2 = [3.0_pr, 1.0_pr]
      call newton_solve(circle_line, x2, tol=1e-12_pr, max_iter=1, its=its, info=info)
      call assert(info == 1, "newton_solve: info = 1 when the iteration limit is reached")
   end subroutine test_newton_solve

   subroutine test_saturation_F()
      !! At the bubble point the residuals vanish, and the analytical
      !! Jacobian matches central finite differences.
      integer, parameter :: nx = nc + 4
      real(pr) :: X(nx)
      real(pr) :: F(nx)
      real(pr) :: F_noder(nx)
      real(pr) :: dF(nx, nx)
      real(pr) :: Fp(nx)
      real(pr) :: Fm(nx)
      real(pr) :: dx(nx)
      real(pr) :: dPdVz
      real(pr) :: dPdVy
      real(pr) :: h
      integer :: ns
      integer :: j

      X = [log(bub%y/z), log(bub%Vx), log(bub%Vy), log(T), log(bub%P)]
      ns = nc + 3

      call saturation_F(model, z, X, ns, X(ns), F, dF, dPdVz, dPdVy)
      call assert(maxval(abs(F(:nc))) < 1e-6_pr, "saturation_F: isofugacity residuals vanish at the bubble point")
      call assert(abs(F(nc+1)) < 1e-10_pr, "saturation_F: material balance residual vanishes")
      call assert(abs(F(nc+2)) < 1e-4_pr*bub%P, "saturation_F: pressure equality at the bubble point")
      call assert(abs(F(nc+4)) == 0, "saturation_F: specification residual vanishes")
      call assert(dPdVz < 0 .and. dPdVy < 0, "saturation_F: mechanically stable phases (dP/dV < 0)")

      ! Without the Jacobian the same residuals are returned
      call saturation_F(model, z, X, ns, X(ns), F_noder, dPdVz=dPdVz, dPdVy=dPdVy)
      call assert(maxval(abs(F_noder - F)) < 1e-12_pr, "saturation_F: residuals independent of dF request")

      do j = 1, nx
         h = 1e-6_pr*max(abs(X(j)), 1.0_pr)
         dx = 0
         dx(j) = h
         call saturation_F(model, z, X + dx, ns, X(ns), Fp, dPdVz=dPdVz, dPdVy=dPdVy)
         call saturation_F(model, z, X - dx, ns, X(ns), Fm, dPdVz=dPdVz, dPdVy=dPdVy)
         call assert(maxval(abs((Fp - Fm)/(2*h) - dF(:, j))) < 1e-4_pr*max(1.0_pr, maxval(abs(dF(:, j)))), &
            "saturation_F: analytical Jacobian column matches finite differences")
      end do
   end subroutine test_saturation_F

   subroutine test_solve_VxVyT()
      !! Newton on the log-volume system from a perturbed bubble point
      !! recovers the bubble pressure at the specified temperature.
      integer, parameter :: nx = nc + 4
      real(pr) :: X(nx)
      integer :: ns
      integer :: its

      X = [log(bub%y/z) + [0.05_pr, -0.05_pr], log(bub%Vx*1.01_pr), log(bub%Vy*1.1_pr), log(T), log(bub%P*1.05_pr)]
      ns = nc + 3
      call solve_VxVyT(model, z, X, ns, log(T), tol=1e-9_pr, max_iterations=50, its=its)
      call assert(its < 50, "solve_VxVyT converges")
      call assert(abs(exp(X(nc+4)) - bub%P) < 1e-4_pr, "solve_VxVyT recovers the bubble pressure")
      call assert(maxval(abs(z*exp(X(:nc)) - bub%y)) < 1e-6_pr, "solve_VxVyT recovers the incipient vapor")
      call assert(abs(exp(X(nc+3)) - T) < 1e-10_pr, "solve_VxVyT keeps the specified temperature")
   end subroutine test_solve_VxVyT

   subroutine test_saturation_TP_stable_roots()
      !! At the bubble point the stable roots are the liquid one for the feed
      !! and the vapor one for the incipient phase, so the stable-root
      !! residuals (any unknown `kind`) equal the bubble residuals.
      real(pr) :: X(nc+2)
      real(pr) :: F_bubble(nc+2)
      real(pr) :: F_stable(nc+2)
      real(pr) :: dF_bubble(nc+2, nc+2)
      real(pr) :: dF_stable(nc+2, nc+2)
      real(pr) :: dFdS(nc+2)

      X = [log(bub%y/z), log(T), log(bub%P)]
      call saturation_TP(model, "bubble", z, X, nc+1, X(nc+1), F_bubble, dF_bubble, dFdS)
      call saturation_TP(model, "stable", z, X, nc+1, X(nc+1), F_stable, dF_stable, dFdS)
      call assert(maxval(abs(F_bubble(:nc))) < 1e-6_pr, "saturation_TP: bubble residuals vanish")
      call assert(maxval(abs(F_stable - F_bubble)) < 1e-10_pr, "saturation_TP: stable roots reproduce the bubble residuals")
      call assert(maxval(abs(dF_stable - dF_bubble)) < 1e-8_pr, "saturation_TP: stable roots reproduce the bubble Jacobian")
      call assert(dFdS(nc+2) == -1 .and. all(dFdS(:nc+1) == 0), "saturation_TP: dF/dS selects the specification row")
   end subroutine test_saturation_TP_stable_roots
end program test_saturation_point_solvers
