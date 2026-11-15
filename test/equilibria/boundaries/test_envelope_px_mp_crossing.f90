program test_envelope_px_mp_crossing
   !! Two-phase lines traced with the multiphase Px solver `px_envelope`
   !! (one main phase, one incipient phase), seeded from converged bubble
   !! points.
   !!
   !! - Ternary: the line crosses the critical point. The recorded critical
   !!   point matches an independent critical-point calculation, the phase
   !!   kinds are swapped after it and every point satisfies isofugacity
   !!   and the material balance. Default ("stable") kinds give the same CP.
   !! - Binary towards the critical point: the line stops at it.
   !! - Binary with `max_pressure`: no point above the limit.
   !! - Binary towards the pure heavy component: the step is limited so
   !!   that \(z\) stays in [0, 1] and the line ends at the pure-component
   !!   vapor pressure.
   !! - Formatted output (`write`) has one row per point.
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, &
      saturation_pressure, critical_point, spec_CP, PXEnvelMP, px_envelope
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Multiphase Px solver: two-phase lines and critical points")

   call test_ternary_crossing()
   call test_binary()

contains

   subroutine test_ternary_crossing()
      integer, parameter :: nc = 3
      real(pr), parameter :: T = 350._pr
      type(CubicEoS) :: model
      type(EquilibriumState) :: sat
      type(EquilibriumState) :: cp
      type(PXEnvelMP) :: env
      type(PXEnvelMP) :: env_stable
      real(pr) :: z0(nc)
      real(pr) :: zi(nc)
      real(pr) :: z(nc)
      real(pr) :: x_l0(1, nc)
      real(pr) :: a
      real(pr) :: lnphi_x(nc)
      real(pr) :: lnphi_w(nc)
      character(len=14) :: kinds_x(1)
      character(len=14) :: kind_w
      integer :: i
      integer :: n
      integer :: n_before
      logical :: ok

      model = PengRobinson76( &
         [190.564_pr, 369.83_pr, 507.6_pr], &
         [45.99_pr, 42.48_pr, 30.25_pr], &
         [0.0115478_pr, 0.1523_pr, 0.3013_pr])
      z0 = [0._pr, 0.5_pr, 0.5_pr]
      zi = [1._pr, 0._pr, 0._pr]

      a = 0.2_pr
      z = a*zi + (1 - a)*z0
      sat = saturation_pressure(model, z, T=T, kind="bubble", P0=50._pr)
      call assert(sat%kind == "bubble", "Ternary Px MP: seeding bubble point converged")

      x_l0(1, :) = z
      kinds_x = "liquid"
      kind_w = "vapor"

      ! X = [lnK(3), beta, lnP, alpha]: specify alpha (position 6)
      env = px_envelope(model, z0, zi, 1, T, x_l0, sat%y, [1._pr], sat%P, a, &
         nc + 3, 0.01_pr, 0._pr, points=80, kinds_x=kinds_x, kind_w=kind_w)
      n = size(env%points)

      call assert(n == 80, "Ternary Px MP: all requested points computed")
      call assert(abs(env%points(1)%P - sat%P) < 1e-6_pr*sat%P, "Ternary Px MP: first point is the seed")

      ! Critical point
      call assert(size(env%Pc) == 1 .and. size(env%ac) == 1, "Ternary Px MP: one critical point recorded")
      cp = critical_point(model, z0, zi, spec=spec_CP%T, S=log(T), max_iters=1000, &
         a0=env%ac(1), P0=env%Pc(1), T0=T)
      call assert(abs(cp%x(1) - env%ac(1)) < 2e-3_pr, "Ternary Px MP: critical alpha matches an independent CP")
      call assert(abs(cp%P - env%Pc(1)) < 0.2_pr, "Ternary Px MP: critical pressure matches an independent CP")

      ! Kinds swapped after the critical point
      n_before = count(env%alpha < env%ac(1))
      ok = all([(env%points(i)%kinds_x(1) == "liquid" .and. env%points(i)%kind_w == "vapor", i=1, n_before)])
      call assert(ok, "Ternary Px MP: liquid main phase before the CP")
      ok = all([(env%points(i)%kinds_x(1) == "vapor" .and. env%points(i)%kind_w == "liquid", i=n_before+1, n)])
      call assert(ok .and. n_before < n, "Ternary Px MP: kinds swapped after the CP")

      do i = 1, n
         z = env%alpha(i)*zi + (1 - env%alpha(i))*z0
         call assert(maxval(abs(env%points(i)%betas(1)*env%points(i)%x_l(1, :) - z)) < 1e-10_pr, &
            "Ternary Px MP: material balance with beta=1")
         call assert(abs(env%points(i)%T - T) < 1e-12_pr, "Ternary Px MP: temperature fixed")
         call model%lnphi_pt(env%points(i)%x_l(1, :), env%points(i)%P, T, &
            root_type=env%points(i)%kinds_x(1), lnPhi=lnphi_x)
         call model%lnphi_pt(env%points(i)%w, env%points(i)%P, T, &
            root_type=env%points(i)%kind_w, lnPhi=lnphi_w)
         call assert(maxval(abs(log(env%points(i)%x_l(1, :)) + lnphi_x &
            - log(env%points(i)%w) - lnphi_w)) < 1e-7_pr, &
            "Ternary Px MP: isofugacity at every point")
         if (i > 1) call assert(env%alpha(i) > env%alpha(i-1), "Ternary Px MP: alpha increases monotonically")
      end do

      ! Without kinds the "stable" roots are used: same critical point
      env_stable = px_envelope(model, z0, zi, 1, T, x_l0, sat%y, [1._pr], sat%P, a, &
         nc + 3, 0.01_pr, 0._pr, points=80)
      call assert(size(env_stable%Pc) == 1, "Ternary Px MP (stable kinds): critical point recorded")
      call assert(abs(env_stable%Pc(1) - env%Pc(1)) < 0.05_pr .and. abs(env_stable%ac(1) - env%ac(1)) < 1e-3_pr, &
         "Ternary Px MP (stable kinds): same critical point")
      call assert(all(env_stable%points%kind_w == "stable"), "Ternary Px MP (stable kinds): stable root kind kept")
   end subroutine test_ternary_crossing

   subroutine test_binary()
      integer, parameter :: nc = 2
      real(pr), parameter :: T = 300._pr
      type(CubicEoS) :: model
      type(EquilibriumState) :: sat
      type(EquilibriumState) :: cp
      type(EquilibriumState) :: psat
      type(PXEnvelMP) :: env
      real(pr) :: z0(nc)
      real(pr) :: zi(nc)
      real(pr) :: z(nc)
      real(pr) :: x_l0(1, nc)
      real(pr) :: a
      character(len=14) :: kinds_x(1)
      character(len=14) :: kind_w
      integer :: n
      integer :: i
      integer :: u
      integer :: ios
      integer :: n_rows
      real(pr) :: row(3)
      character(len=1024) :: line

      ! Methane injected into (almost pure) propane
      model = PengRobinson76([190.564_pr, 369.83_pr], [45.99_pr, 42.48_pr], [0.0115478_pr, 0.1523_pr])
      z0 = [1e-3_pr, 1 - 1e-3_pr]
      zi = [1._pr, 0._pr]

      a = 0.05_pr
      z = a*zi + (1 - a)*z0
      sat = saturation_pressure(model, z, T=T, kind="bubble", P0=20._pr)
      x_l0(1, :) = z
      kinds_x = "liquid"
      kind_w = "vapor"

      ! ------------------------------------------------------------------------
      ! Towards the binary critical point: stops at it
      ! ------------------------------------------------------------------------
      env = px_envelope(model, z0, zi, 1, T, x_l0, sat%y, [1._pr], sat%P, a, &
         nc + 3, 0.01_pr, 0._pr, points=500, kinds_x=kinds_x, kind_w=kind_w)
      n = size(env%points)
      call assert(n > 5 .and. n < 500, "Binary Px MP: line stops before the maximum number of points")
      cp = critical_point(model, z0, zi, spec=spec_CP%T, S=log(T), max_iters=1000, &
         a0=env%alpha(n), P0=env%points(n)%P, T0=T)
      call assert(abs(cp%T - T) < 1e-6_pr, "Binary Px MP: independent CP converged")
      call assert(abs(env%points(n)%P - cp%P) < 1._pr, "Binary Px MP: line ends at the binary critical pressure")
      call assert(all(env%points%P <= cp%P + 0.1_pr), "Binary Px MP: no point beyond the critical pressure")

      ! Formatted output: one row per point, starting with alpha, P, T
      open(newunit=u, status="scratch", form="formatted", action="readwrite")
      call env%write(u)
      rewind(u)
      n_rows = 0
      do
         read(u, "(A)", iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         n_rows = n_rows + 1
         if (n_rows == 1) read(line, *) row
      end do
      close(u)
      call assert(n_rows == n, "PXEnvelMP write: one row per point")
      call assert(abs(row(1) - env%alpha(1)) < 1e-4_pr*abs(env%alpha(1)) &
         .and. abs(row(2) - env%points(1)%P) < 1e-4_pr*env%points(1)%P &
         .and. abs(row(3) - T) < 1e-4_pr*T, "PXEnvelMP write: first row is alpha, P, T")

      ! ------------------------------------------------------------------------
      ! Pressure limit
      ! ------------------------------------------------------------------------
      env = px_envelope(model, z0, zi, 1, T, x_l0, sat%y, [1._pr], sat%P, a, &
         nc + 3, 0.01_pr, 0._pr, points=500, kinds_x=kinds_x, kind_w=kind_w, max_pressure=50._pr)
      n = size(env%points)
      call assert(all(env%points%P <= 50._pr), "Binary Px MP: max_pressure respected")
      call assert(env%points(n)%P > 40._pr, "Binary Px MP: line traced up to the pressure limit")

      ! ------------------------------------------------------------------------
      ! Towards pure propane: z is kept inside [0, 1]
      ! ------------------------------------------------------------------------
      env = px_envelope(model, z0, zi, 1, T, x_l0, sat%y, [1._pr], sat%P, a, &
         nc + 3, -0.01_pr, 0._pr, points=500, kinds_x=kinds_x, kind_w=kind_w)
      n = size(env%points)
      do i = 1, n
         z = env%alpha(i)*zi + (1 - env%alpha(i))*z0
         call assert(all(z >= -1e-12_pr) .and. all(z <= 1 + 1e-12_pr), "Binary Px MP: z stays inside [0, 1]")
      end do
      psat = saturation_pressure(model, [0._pr, 1._pr], T=T, kind="bubble", P0=10._pr)
      call assert(abs(env%points(n)%P - psat%P) < 1e-5_pr*psat%P, &
         "Binary Px MP: line ends at the pure propane vapor pressure")
   end subroutine test_binary
end program test_envelope_px_mp_crossing
