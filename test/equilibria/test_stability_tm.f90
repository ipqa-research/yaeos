program test_stability_tm
   !! Tangent plane distance (`tm`) and its minimization (`min_tpd`) for
   !! residual (ArModel) and excess-Gibbs (GeModel) models:
   !! tm(w = z) = 0, analytical gradient vs finite differences, sign of the
   !! minimum for stable/unstable feeds, consistency with flash, and the
   !! `all_minima` output. Also checks the formatted output of an
   !! `EquilibriumState`.
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, flash, &
      UNIFAC, Groups, setup_unifac
   use yaeos__equilibria_stability, only: tm, min_tpd
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   class(ArModel), allocatable :: model
   real(pr) :: z(nc)
   real(pr) :: T

   write(*, *) test_title("Tangent plane distance")

   ! Methane / n-butane
   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]
   T = 300.0_pr

   call test_ar_tm()
   call test_ar_min_tpd()
   call test_ge()
   call test_state_output()

contains

   subroutine test_ar_tm()
      real(pr) :: w(nc)
      real(pr) :: d(nc)
      real(pr) :: lnphi_z(nc)
      real(pr) :: grad(nc)
      real(pr) :: dw(nc)
      real(pr) :: num
      real(pr) :: h
      real(pr) :: V
      real(pr) :: P
      integer :: i

      P = 40.0_pr
      call assert(abs(tm(model, z, z, P, T)) < 1e-12_pr, "ArModel: tm(w = z) = 0")

      ! Gradient vs central differences, with the feed d_i given
      call model%lnphi_pt(z, P, T, V=V, root_type="stable", lnPhi=lnphi_z)
      d = log(z) + lnphi_z
      w = [0.9_pr, 0.1_pr]
      call assert(abs(tm(model, z, w, P, T, d=d) - tm(model, z, w, P, T)) < 1e-12_pr, &
         "ArModel: tm with the feed d_i given equals tm computing them")
      num = tm(model, z, w, P, T, d=d, dtpd=grad)
      do i = 1, nc
         h = 1e-6_pr
         dw = 0
         dw(i) = h
         num = (tm(model, z, w + dw, P, T, d=d) - tm(model, z, w - dw, P, T, d=d))/(2*h)
         call assert(abs(num - grad(i)) < 1e-6_pr, "ArModel: tm gradient matches finite differences")
      end do
   end subroutine test_ar_tm

   subroutine test_ar_min_tpd()
      type(EquilibriumState) :: fr
      real(pr) :: w(nc)
      real(pr) :: mintpd
      real(pr) :: all_minima(nc, nc + 1)
      integer :: i

      ! Unstable feed in the two-phase region
      call min_tpd(model, z, 40.0_pr, T, mintpd, w, all_minima=all_minima)
      call assert(mintpd < -1e-3_pr, "40 bar: feed is unstable (negative tm minimum)")
      call assert(abs(sum(w) - 1) < 1e-12_pr, "min_tpd: trial composition is normalized")
      call assert(abs(tm(model, z, w, 40.0_pr, T) - mintpd) < 1e-10_pr, "min_tpd: returned value is tm at w")
      call assert(abs(minval(all_minima(:, nc + 1)) - mintpd) < 1e-14_pr, "all_minima holds the global minimum")
      do i = 1, nc
         call assert(abs(sum(all_minima(i, :nc)) - 1) < 1e-12_pr, "all_minima: every trial phase is normalized")
         call assert(abs(tm(model, z, all_minima(i, :nc), 40.0_pr, T) - all_minima(i, nc + 1)) < 1e-10_pr, &
            "all_minima: stored values are tm at the stored compositions")
      end do
      ! The most negative trial phase is the vapor-like (methane-rich) one;
      ! it lies close to the flash vapor.
      fr = flash(model, z, T, p_spec=40.0_pr)
      call assert(w(1) > z(1), "40 bar: unstable trial phase is methane-rich")
      call assert(abs(w(1) - fr%y(1)) < 0.1_pr, "40 bar: trial phase close to the flash vapor")

      ! Stable compressed liquid and stable vapor (below the dew pressure)
      call min_tpd(model, z, 200.0_pr, T, mintpd, w)
      call assert(mintpd > -1e-8_pr, "200 bar: compressed liquid is stable")
      call min_tpd(model, z, 3.0_pr, T, mintpd, w)
      call assert(mintpd > -1e-8_pr, "3 bar: vapor below the dew pressure is stable")
   end subroutine test_ar_min_tpd

   subroutine test_ge()
      !! UNIFAC liquid-liquid system of `test_flash_ge`.
      integer, parameter :: ncg = 3
      type(UNIFAC) :: ge
      type(Groups) :: molecules(ncg)
      real(pr) :: n(ncg)
      real(pr) :: w(ncg)
      real(pr) :: grad(ncg)
      real(pr) :: dw(ncg)
      real(pr) :: all_minima(ncg, ncg + 1)
      real(pr) :: mintpd
      real(pr) :: val
      real(pr) :: num
      real(pr) :: h
      integer :: i

      n = [0.2_pr, 0.7_pr, 0.1_pr]
      molecules(1)%groups_ids = [1, 42]
      molecules(1)%number_of_groups = [1, 1]
      molecules(2)%groups_ids = [1, 2]
      molecules(2)%number_of_groups = [2, 4]
      molecules(3)%groups_ids = [16]
      molecules(3)%number_of_groups = [1]
      ge = setup_unifac(molecules)

      call assert(abs(tm(ge, n, n, 1.0_pr, 250.0_pr)) < 1e-12_pr, "GeModel: tm(w = z) = 0")

      call min_tpd(ge, n, 1.0_pr, 250.0_pr, mintpd, w, all_minima=all_minima)
      call assert(mintpd < 0, "GeModel: feed is unstable")
      call assert(abs(tm(ge, n, w, 1.0_pr, 250.0_pr) - mintpd) < 1e-10_pr, "GeModel: min_tpd value is tm at w")
      call assert(abs(minval(all_minima(:, ncg + 1)) - mintpd) < 1e-14_pr, "GeModel: all_minima holds the minimum")

      w = [0.3_pr, 0.3_pr, 0.4_pr]
      val = tm(ge, n, w, 1.0_pr, 250.0_pr, dtpd=grad)
      do i = 1, ncg
         h = 1e-6_pr
         dw = 0
         dw(i) = h
         num = (tm(ge, n, w + dw, 1.0_pr, 250.0_pr) - tm(ge, n, w - dw, 1.0_pr, 250.0_pr))/(2*h)
         call assert(abs(num - grad(i)) < 1e-6_pr, "GeModel: tm gradient matches finite differences")
      end do
   end subroutine test_ge

   subroutine test_state_output()
      !! Formatted (list-directed) output of an EquilibriumState holds the
      !! kind, T, P, beta, x and y, in that order.
      type(EquilibriumState) :: fr
      character(len=500) :: line
      character(len=14) :: kind
      real(pr) :: Tr
      real(pr) :: Pr_
      real(pr) :: beta
      real(pr) :: x(nc)
      real(pr) :: y(nc)

      fr = flash(model, z, T, p_spec=40.0_pr)
      write(line, *) fr
      read(line, *) kind, Tr, Pr_, beta, x, y
      call assert(kind == fr%kind, "EquilibriumState output: kind")
      call assert(abs(Tr - fr%T) < 1e-10_pr .and. abs(Pr_ - fr%P) < 1e-10_pr, "EquilibriumState output: T and P")
      call assert(abs(beta - fr%beta) < 1e-12_pr, "EquilibriumState output: beta")
      call assert(maxval(abs(x - fr%x)) < 1e-12_pr .and. maxval(abs(y - fr%y)) < 1e-12_pr, &
         "EquilibriumState output: phase compositions")
   end subroutine test_state_output
end program test_stability_tm
