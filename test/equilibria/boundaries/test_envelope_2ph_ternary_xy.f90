program test_envelope_2ph_ternary_xy
   !! Two-phase Px and Tx envelopes of a ternary system (methane added to a
   !! propane + n-hexane mixture) started from a dew point.
   !!
   !! The lines start on the dew branch, cross the mixture critical point and
   !! continue on the bubble branch. Checked invariants:
   !!
   !! - The first point reproduces the seeding saturation point.
   !! - Every point satisfies isofugacity and the fixed variable (T or P).
   !! - Points agree with independently computed saturation points.
   !! - The detected critical point agrees with an independent critical-point
   !!   calculation at the fixed variable, and the point kind switches there.
   !! - `stop_conditions` ends the line where requested.
   !! - A `"liquid-liquid"` first point (stable roots) traces the same line.
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, &
      saturation_pressure, saturation_temperature, critical_point, spec_CP, &
      PXEnvel2, px_envelope_2ph, TXEnvel2, tx_envelope_2ph
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 3
   type(CubicEoS) :: model
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)

   write(*, *) test_title("Two-phase Px/Tx envelopes of a ternary from a dew point")

   model = PengRobinson76( &
      [190.564_pr, 369.83_pr, 507.6_pr], &
      [45.99_pr, 42.48_pr, 30.25_pr], &
      [0.0115478_pr, 0.1523_pr, 0.3013_pr])
   z0 = [0._pr, 0.5_pr, 0.5_pr]
   zi = [1._pr, 0._pr, 0._pr]

   call test_px()
   call test_tx()
   call test_liquid_liquid_kind()

contains

   real(pr) function isofugacity_error(x, y, P, T, kind_x, kind_y)
      !! Maximum absolute residual of \(\ln x_i\hat\phi_i^x = \ln y_i\hat\phi_i^y\).
      real(pr), intent(in) :: x(nc)
      real(pr), intent(in) :: y(nc)
      real(pr), intent(in) :: P
      real(pr), intent(in) :: T
      character(len=*), intent(in) :: kind_x
      character(len=*), intent(in) :: kind_y
      real(pr) :: lnphi_x(nc)
      real(pr) :: lnphi_y(nc)
      call model%lnphi_pt(x, P, T, root_type=kind_x, lnPhi=lnphi_x)
      call model%lnphi_pt(y, P, T, root_type=kind_y, lnPhi=lnphi_y)
      isofugacity_error = maxval(abs(log(x) + lnphi_x - log(y) - lnphi_y), mask=x > 0 .and. y > 0)
   end function isofugacity_error

   logical function stop_px(X, ns, S, dS, dXdS, iterations)
      !! Stop the Px line once the injected fraction leaves [0.02, 1].
      real(pr), intent(in out) :: X(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: dS
      real(pr), intent(in out) :: dXdS(:)
      integer, intent(in) :: iterations
      stop_px = X(nc+2) + dXdS(nc+2)*dS < 0.02_pr .or. X(nc+2) + dXdS(nc+2)*dS > 1
   end function stop_px

   subroutine test_px()
      real(pr), parameter :: T = 350._pr
      type(EquilibriumState) :: sat
      type(EquilibriumState) :: check
      type(EquilibriumState) :: cp
      type(PXEnvel2) :: env
      real(pr) :: a
      real(pr) :: z(nc)
      real(pr) :: err
      integer :: i
      integer :: n
      integer :: i_bub
      integer :: n_dew
      integer :: n_bubble
      logical :: kinds_ok

      a = 0.9_pr
      z = a*zi + (1 - a)*z0
      sat = saturation_pressure(model, z, T=T, kind="dew", P0=10._pr)
      call assert(sat%kind == "dew", "Px: seeding dew point converged")

      env = px_envelope_2ph(model, z0, a, zi, sat, delta_0=0.01_pr, points=200, &
         stop_conditions=stop_px)
      n = size(env%points)

      call assert(n > 10 .and. n < 200, "Px: line traced and stopped by stop_conditions")
      call assert(abs(env%points(1)%P - sat%P) < 1e-4_pr*sat%P, "Px: first point is the seeding dew point")
      call assert(abs(env%alpha(1) - a) < 1e-12_pr, "Px: first point at the seeding alpha")
      call assert(env%alpha(n) >= 0.0_pr .and. env%alpha(n) < 0.15_pr, &
         "Px: stop_conditions ends the line near alpha = 0")

      ! Kinds: dew points, then bubble points after the critical point
      n_dew = count(env%points%kind == "dew")
      n_bubble = count(env%points%kind == "bubble")
      call assert(n_dew > 0 .and. n_bubble > 0 .and. n_dew + n_bubble == n, &
         "Px: line has a dew and a bubble branch")
      kinds_ok = all(env%points(:n_dew)%kind == "dew")
      call assert(kinds_ok, "Px: all dew points come before the bubble points")

      ! Isofugacity and fixed temperature on every point
      do i = 1, n
         call assert(abs(env%points(i)%T - T) < 1e-10_pr, "Px: temperature is fixed along the line")
         err = isofugacity_error(env%points(i)%x, env%points(i)%y, env%points(i)%P, T, "liquid", "vapor")
         call assert(err < 1e-4_pr, "Px: isofugacity at every envelope point")
      end do

      ! Critical point: one, between the branches, matching an independent CP
      call assert(size(env%cps) == 1, "Px: one critical point found")
      cp = critical_point(model, z0, zi, spec=spec_CP%T, S=log(T), max_iters=1000, &
         a0=env%cps(1)%alpha, P0=env%cps(1)%P, T0=T)
      call assert(cp%x(1) < env%alpha(n_dew) .and. cp%x(1) > env%alpha(n_dew + 1), &
         "Px: the true critical alpha lies between the last dew and first bubble point")
      call assert(abs(cp%T - T) < 1e-6_pr, "Px: independent CP converged at the envelope temperature")
      call assert(abs(cp%x(1) - env%cps(1)%alpha) < 0.02_pr, "Px: critical alpha agrees with an independent CP")
      call assert(abs(cp%P - env%cps(1)%P) < 2._pr, "Px: critical pressure agrees with an independent CP")

      ! A bubble point well inside the bubble branch matches saturation_pressure
      i_bub = n_dew + (n_bubble/2)
      check = saturation_pressure(model, env%points(i_bub)%x, T=T, kind="bubble", P0=env%points(i_bub)%P)
      call assert(abs(check%P - env%points(i_bub)%P) < 1e-3_pr*check%P, &
         "Px: bubble point on the line matches saturation_pressure")
      ! A dew point (incipient phase x, main phase y = z)
      check = saturation_pressure(model, env%points(2)%y, T=T, kind="dew", P0=env%points(2)%P)
      call assert(abs(check%P - env%points(2)%P) < 1e-3_pr*check%P, &
         "Px: dew point on the line matches saturation_pressure")
   end subroutine test_px

   logical function stop_tx(X, ns, S, dS, dXdS, iterations)
      !! Stop the Tx line once the injected fraction leaves [0.02, 0.98].
      real(pr), intent(in out) :: X(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: S
      real(pr), intent(in out) :: dS
      real(pr), intent(in out) :: dXdS(:)
      integer, intent(in) :: iterations
      stop_tx = X(nc+2) + dXdS(nc+2)*dS < 0.02_pr .or. X(nc+2) + dXdS(nc+2)*dS > 0.98_pr
   end function stop_tx

   subroutine test_tx()
      real(pr), parameter :: P = 80._pr
      type(EquilibriumState) :: sat
      type(EquilibriumState) :: check
      type(EquilibriumState) :: cp
      type(TXEnvel2) :: env
      real(pr) :: a
      real(pr) :: z(nc)
      real(pr) :: err
      integer :: i
      integer :: n
      integer :: i_bub
      integer :: n_dew
      integer :: n_bubble

      a = 0.4_pr
      z = a*zi + (1 - a)*z0
      sat = saturation_temperature(model, z, P=P, kind="dew", T0=450._pr)
      call assert(sat%kind == "dew", "Tx: seeding dew point converged")

      env = tx_envelope_2ph(model, z0, a, zi, sat, delta_0=-0.01_pr, points=300, &
         stop_conditions=stop_tx)
      n = size(env%points)

      call assert(n > 10 .and. n < 300, "Tx: line traced and stopped by stop_conditions")
      call assert(abs(env%points(1)%T - sat%T) < 0.01_pr, "Tx: first point is the seeding dew point")
      call assert(env%alpha(n) > 0.0_pr .and. env%alpha(n) < 1.0_pr, "Tx: line ends inside 0 < alpha < 1")

      n_dew = count(env%points%kind == "dew")
      n_bubble = count(env%points%kind == "bubble")
      call assert(n_dew > 0 .and. n_bubble > 0 .and. n_dew + n_bubble == n, &
         "Tx: line has a dew and a bubble branch")
      call assert(all(env%points(:n_dew)%kind == "dew"), "Tx: all dew points come before the bubble points")

      do i = 1, n
         call assert(abs(env%points(i)%P - P) < 1e-10_pr, "Tx: pressure is fixed along the line")
         err = isofugacity_error(env%points(i)%x, env%points(i)%y, P, env%points(i)%T, "liquid", "vapor")
         call assert(err < 1e-4_pr, "Tx: isofugacity at every envelope point")
         if (i > 1) call assert(abs(env%points(i)%T - env%points(i-1)%T) < 25._pr, &
            "Tx: temperature is continuous between consecutive points")
      end do

      ! Only the first critical point is checked: further spurious CPs with
      ! alpha = 0 can be reported on the bubble branch (see the interpolation
      ! on X(ns) in the Tx detect_critical).
      call assert(size(env%cps) >= 1, "Tx: critical point found")
      cp = critical_point(model, z0, zi, spec=spec_CP%P, S=log(P), max_iters=1000, &
         a0=env%cps(1)%alpha, T0=env%cps(1)%T)
      call assert(cp%x(1) < env%alpha(n_dew) .and. cp%x(1) > env%alpha(n_dew + 1), &
         "Tx: the true critical alpha lies between the last dew and first bubble point")
      call assert(abs(cp%P - P) < 1e-6_pr, "Tx: independent CP converged at the envelope pressure")
      call assert(abs(cp%x(1) - env%cps(1)%alpha) < 0.02_pr, "Tx: critical alpha agrees with an independent CP")
      call assert(abs(cp%T - env%cps(1)%T) < 2._pr, "Tx: critical temperature agrees with an independent CP")

      i_bub = n_dew + (n_bubble/2)
      check = saturation_temperature(model, env%points(i_bub)%x, P=P, kind="bubble", T0=env%points(i_bub)%T)
      call assert(abs(check%T - env%points(i_bub)%T) < 0.05_pr, &
         "Tx: bubble point on the line matches saturation_temperature")
   end subroutine test_tx

   subroutine test_liquid_liquid_kind()
      !! A first point labelled "liquid-liquid" uses the stable volume roots.
      !! Started from a VLE bubble point it must trace the same line as the
      !! bubble-labelled start, with points labelled "saturation".
      real(pr), parameter :: T = 350._pr
      real(pr), parameter :: P = 80._pr
      type(EquilibriumState) :: sat
      type(EquilibriumState) :: ll
      type(PXEnvel2) :: px_bub
      type(PXEnvel2) :: px_ll
      type(TXEnvel2) :: tx_bub
      type(TXEnvel2) :: tx_ll
      real(pr) :: a
      real(pr) :: z(nc)
      integer :: i

      a = 0.2_pr
      z = a*zi + (1 - a)*z0

      sat = saturation_pressure(model, z, T=T, kind="bubble", P0=50._pr)
      ll = sat
      ll%kind = "liquid-liquid"
      px_bub = px_envelope_2ph(model, z0, a, zi, sat, delta_0=0.01_pr, points=5)
      px_ll = px_envelope_2ph(model, z0, a, zi, ll, delta_0=0.01_pr, points=5)
      call assert(size(px_ll%points) == 5, "Px liquid-liquid start: line traced")
      call assert(all(px_ll%points%kind == "saturation"), "Px liquid-liquid start: points labelled saturation")
      do i = 1, 5
         call assert(abs(px_ll%points(i)%P - px_bub%points(i)%P) < 1e-6_pr*px_bub%points(i)%P, &
            "Px liquid-liquid start: same line as the bubble start")
      end do

      a = 0.53_pr
      z = a*zi + (1 - a)*z0
      sat = saturation_temperature(model, z, P=P, kind="bubble", T0=300._pr)
      ll = sat
      ll%kind = "liquid-liquid"
      tx_bub = tx_envelope_2ph(model, z0, a, zi, sat, delta_0=0.01_pr, points=5)
      tx_ll = tx_envelope_2ph(model, z0, a, zi, ll, delta_0=0.01_pr, points=5)
      call assert(size(tx_ll%points) == 5, "Tx liquid-liquid start: line traced")
      call assert(all(tx_ll%points%kind == "saturation"), "Tx liquid-liquid start: points labelled saturation")
      do i = 1, 5
         call assert(abs(tx_ll%points(i)%T - tx_bub%points(i)%T) < 1e-6_pr*tx_bub%points(i)%T, &
            "Tx liquid-liquid start: same line as the bubble start")
      end do
   end subroutine test_liquid_liquid_kind
end program test_envelope_2ph_ternary_xy
