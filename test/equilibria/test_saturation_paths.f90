program test_saturation_paths
   !! Bubble, dew and liquid-liquid saturation points through every route of
   !! `saturation_pressure` and `saturation_temperature`: default and user
   !! initial guesses, the full-Newton fallback (forced through the public
   !! `iters_first_step` module variable) and failure detection.
   !! Converged points are checked for isofugacity and normalization, against
   !! flash calculations, and `saturation_temperature` is checked to be the
   !! inverse of `saturation_pressure`.
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, flash, &
      saturation_pressure, saturation_temperature
   use yaeos__equilibria_saturation_points, only: iters_first_step
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: bub
   type(EquilibriumState) :: dew
   type(EquilibriumState) :: res
   type(EquilibriumState) :: fr
   real(pr) :: z(nc)
   real(pr) :: T
   integer :: default_first_step

   write(*, *) test_title("Saturation point routes")

   ! Methane / n-butane
   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]
   T = 300.0_pr
   default_first_step = iters_first_step

   ! ===========================================================================
   ! Bubble pressure (default initialization)
   ! ---------------------------------------------------------------------------
   bub = saturation_pressure(model, z, T, kind="bubble")
   call assert(bub%kind == "bubble", "bubble pressure converges")
   call check_point(bub, "bubble pressure", 1e-6_pr)
   call assert(bub%beta == 0, "bubble point has beta = 0")
   call assert(maxval(abs(bub%x - z)) < 1e-12_pr, "bubble point: liquid is the feed")
   call assert(abs(bub%P - 57.3726_pr) < 1e-3_pr, "bubble pressure of 30% methane at 300 K")

   ! Flash just below / above the bubble pressure
   fr = flash(model, z, T, p_spec=bub%P*(1 - 1e-3_pr))
   call assert(fr%kind == "split", "flash slightly below the bubble pressure splits")
   call assert(fr%beta < 1e-2_pr, "flash slightly below the bubble pressure: little vapor")
   call assert(maxval(abs(fr%y - bub%y)) < 1e-2_pr, "flash vapor approaches the incipient bubble")
   fr = flash(model, z, T, p_spec=bub%P*(1 + 1e-3_pr))
   call assert(fr%kind == "failed", "flash slightly above the bubble pressure is single phase")

   ! User initial guess (converged incipient phase): immediate convergence
   res = saturation_pressure(model, z, T, kind="bubble", p0=bub%P*1.01_pr, y0=bub%y)
   call assert(res%kind == "bubble", "bubble pressure with p0, y0 converges")
   call assert(abs(res%P - bub%P) < 1e-4_pr, "bubble pressure with p0, y0: same pressure")
   call assert(res%iters < 10, "bubble pressure with a good y0 needs few iterations")

   ! Inverse: bubble temperature at the bubble pressure
   res = saturation_temperature(model, z, bub%P, kind="bubble", t0=280.0_pr)
   call assert(res%kind == "bubble", "bubble temperature converges")
   call check_point(res, "bubble temperature", 1e-3_pr)
   call assert(abs(res%T - T) < 1e-3_pr, "bubble temperature inverts the bubble pressure")
   call assert(maxval(abs(res%y - bub%y)) < 1e-3_pr, "bubble temperature: same incipient vapor")
   ! The default saturation_temperature also overwrites the first-step limit,
   ! so restore it for the following checks.
   iters_first_step = default_first_step

   res = saturation_temperature(model, z, bub%P, kind="bubble", t0=290.0_pr, y0=bub%y)
   call assert(abs(res%T - T) < 1e-3_pr, "bubble temperature with y0 inverts the bubble pressure")
   iters_first_step = default_first_step

   ! ===========================================================================
   ! Dew pressure (initial pressure given: see report on P_wilson)
   ! ---------------------------------------------------------------------------
   dew = saturation_pressure(model, z, T, kind="dew", p0=3.0_pr)
   call assert(dew%kind == "dew", "dew pressure converges")
   call check_point(dew, "dew pressure", 1e-3_pr)
   call assert(dew%beta == 1, "dew point has beta = 1")
   call assert(maxval(abs(dew%y - z)) < 1e-12_pr, "dew point: vapor is the feed")
   call assert(dew%P < bub%P, "dew pressure is below the bubble pressure")

   fr = flash(model, z, T, p_spec=dew%P*(1 + 1e-3_pr))
   call assert(fr%kind == "split", "flash slightly above the dew pressure splits")
   call assert(fr%beta > 0.99_pr, "flash slightly above the dew pressure: little liquid")

   res = saturation_temperature(model, z, dew%P, kind="dew", t0=310.0_pr)
   call assert(res%kind == "dew", "dew temperature converges")
   call check_point(res, "dew temperature", 1e-3_pr)
   call assert(abs(res%T - T) < 1e-3_pr, "dew temperature inverts the dew pressure")
   iters_first_step = default_first_step

   ! ===========================================================================
   ! Full-Newton fallback: one scalar step only, then solve_TP
   ! ---------------------------------------------------------------------------
   iters_first_step = 1
   res = saturation_pressure(model, z, T, kind="bubble", p0=50.0_pr)
   call assert(res%kind == "bubble", "fallback bubble pressure converges")
   call check_point(res, "fallback bubble pressure", 1e-6_pr)
   call assert(abs(res%P - bub%P) < 1e-4_pr, "fallback bubble pressure equals the scalar route")
   call assert(maxval(abs(res%y - bub%y)) < 1e-6_pr, "fallback bubble: same incipient vapor")
   call assert(abs(res%Vx - bub%Vx) < 1e-6_pr*bub%Vx, "fallback bubble: same liquid molar volume")

   res = saturation_pressure(model, z, T, kind="dew", p0=3.5_pr)
   call assert(res%kind == "dew", "fallback dew pressure converges")
   call check_point(res, "fallback dew pressure", 1e-6_pr)
   call assert(abs(res%P - dew%P) < 1e-4_pr, "fallback dew pressure equals the scalar route")
   iters_first_step = default_first_step

   ! ===========================================================================
   ! No saturation point: trivial solution flagged as failed
   ! ---------------------------------------------------------------------------
   ! Above the critical temperature of both components, starting from y = z
   res = saturation_pressure(model, z, 500.0_pr, kind="bubble", y0=z)
   call assert(res%kind == "failed", "bubble pressure above both Tc is flagged as failed")
   res = saturation_temperature(model, z, 500.0_pr, kind="bubble", t0=600.0_pr, y0=z)
   call assert(res%kind == "failed", "bubble temperature with a trivial initial guess is flagged as failed")
   iters_first_step = default_first_step

   ! Scalar step not converging from 100 K at 20 bar: the full-Newton
   ! fallback is used and its failure is reported.
   res = saturation_temperature(model, z, 20.0_pr, kind="dew", t0=100.0_pr)
   call assert(res%kind == "failed", "non-converged dew temperature is flagged as failed")
   iters_first_step = default_first_step

   call test_liquid_liquid()

contains

   subroutine check_point(sat, name, tol)
      !! Normalization and isofugacity of a saturation point, with the phases
      !! evaluated at their stored roots. `tol` is the isofugacity tolerance
      !! on \(\ln f_i\).
      type(EquilibriumState), intent(in) :: sat
      character(len=*), intent(in) :: name
      real(pr), intent(in) :: tol
      real(pr) :: lnphi_x(nc)
      real(pr) :: lnphi_y(nc)
      real(pr) :: Vx
      real(pr) :: Vy

      call assert(abs(sum(sat%x) - 1) < 1e-5_pr, name // ": x is normalized")
      call assert(abs(sum(sat%y) - 1) < 1e-5_pr, name // ": y is normalized")
      call assert(maxval(abs(sat%x - sat%y)) > 1e-2_pr, name // ": non-trivial point")

      ! Fugacities recomputed at (T, P) with the expected roots: the
      ! liquid fugacity is too sensitive to V for a VT evaluation at the
      ! stored (previous-iterate) volume to be meaningful.
      call model%lnphi_pt(sat%x, sat%P, sat%T, V=Vx, root_type="liquid", lnPhi=lnphi_x)
      call model%lnphi_pt(sat%y, sat%P, sat%T, V=Vy, root_type="vapor", lnPhi=lnphi_y)
      call assert(abs(Vx - sat%Vx) < 1e-3_pr*Vx, name // ": stored liquid volume is the liquid root")
      call assert(abs(Vy - sat%Vy) < 1e-3_pr*Vy, name // ": stored vapor volume is the vapor root")
      call assert(maxval(abs(log(sat%x) + lnphi_x - log(sat%y) - lnphi_y)) < tol, &
         name // ": isofugacity")
   end subroutine check_point

   subroutine test_liquid_liquid()
      !! Liquid-liquid saturation of the heavy phase of a LL flash: both
      !! saturation routines must return the flash conditions and the
      !! flash light phase as the incipient phase.
      class(ArModel), allocatable :: ll
      type(EquilibriumState) :: lle
      type(EquilibriumState) :: sat
      real(pr) :: kij(nc, nc)
      real(pr) :: zll(nc)
      real(pr) :: Pll

      kij = 0
      kij(1, 2) = 0.2_pr
      kij(2, 1) = 0.2_pr
      ll = PengRobinson76([500.0_pr, 520.0_pr], [30.0_pr, 28.0_pr], [0.3_pr, 0.35_pr], kij)
      zll = [0.5_pr, 0.5_pr]
      Pll = 10.0_pr

      lle = flash(ll, zll, T, p_spec=Pll, k0=[10.0_pr, 0.1_pr])
      call assert(lle%kind == "split", "LL flash splits")
      call assert(lle%Vx < 0.2_pr .and. lle%Vy < 0.2_pr, "LL flash: both phases are liquids")

      sat = saturation_pressure(ll, lle%x, T, kind="liquid-liquid", p0=12.0_pr, y0=lle%y)
      call assert(sat%kind == "liquid-liquid", "liquid-liquid saturation pressure converges")
      call assert(abs(sat%P - Pll) < 1e-2_pr, "LL saturation pressure recovers the flash pressure")
      call assert(maxval(abs(sat%y - lle%y)) < 1e-4_pr, "LL saturation pressure: incipient phase = flash light phase")
      call assert(sat%Vy < 0.2_pr, "LL saturation pressure: incipient phase is a liquid")

      sat = saturation_temperature(ll, lle%x, Pll, kind="liquid-liquid", t0=310.0_pr, y0=lle%y)
      call assert(sat%kind == "liquid-liquid", "liquid-liquid saturation temperature converges")
      call assert(abs(sat%T - T) < 1e-3_pr, "LL saturation temperature recovers the flash temperature")
      call assert(maxval(abs(sat%y - lle%y)) < 1e-3_pr, "LL saturation temperature: incipient phase = flash light phase")
      iters_first_step = default_first_step

      iters_first_step = 1
      sat = saturation_pressure(ll, lle%x, T, kind="liquid-liquid", p0=12.0_pr, y0=lle%y)
      call assert(sat%kind == "liquid-liquid", "fallback LL saturation pressure converges")
      call assert(abs(sat%P - Pll) < 1e-2_pr, "fallback LL saturation pressure recovers the flash pressure")
      call assert(maxval(abs(sat%y - lle%y)) < 1e-6_pr, "fallback LL: incipient phase = flash light phase")
      iters_first_step = default_first_step
   end subroutine test_liquid_liquid
end program test_saturation_paths
