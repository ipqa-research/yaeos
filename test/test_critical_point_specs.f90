program test_critical_point_specs
   !! Mixture critical points with every specification variable (\(\alpha\),
   !! \(V\), \(T\), \(P\)) and with/without initial guesses, the critical line
   !! with a pressure limit, and the recovery of pure-component critical
   !! constants from the EoS (`get_critical_constants`).
   !! Results are checked against the criticality conditions, against each
   !! other, and against the critical point detected on a PT phase envelope.
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, &
      critical_point, critical_line, CriticalLine, spec_CP, &
      saturation_temperature, pt_envelope_2ph, PTEnvel2
   use yaeos__equilibria_critical, only: lambda1, F_critical, get_critical_constants
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   type(CubicEoS) :: model
   type(EquilibriumState) :: cp
   type(EquilibriumState) :: other
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)
   real(pr), parameter :: a_spec = 0.7_pr

   write(*, *) test_title("Critical point specifications")

   ! Methane / n-butane
   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z0 = [1.0_pr, 0.0_pr]
   zi = [0.0_pr, 1.0_pr]

   ! ===========================================================================
   ! Specified alpha, default initial guesses
   ! ---------------------------------------------------------------------------
   cp = critical_point(model, z0, zi, spec_CP%a, S=a_spec, max_iters=1000)
   call assert(cp%kind == "critical", "critical point at specified alpha")
   call assert(maxval(abs(cp%x - [1 - a_spec, a_spec])) < 1e-10_pr, "critical point keeps the specified composition")
   call assert(maxval(abs(cp%y - cp%x)) == 0 .and. cp%Vx == cp%Vy, "critical point: identical phases")
   call check_criticality(cp, "alpha spec")
   call compare_with_envelope(cp)

   ! ===========================================================================
   ! Other specifications must reproduce the same point
   ! ---------------------------------------------------------------------------
   ! Temperature, no initial guesses (alpha starts at 0.5)
   other = critical_point(model, z0, zi, spec_CP%T, S=log(cp%T), max_iters=1000)
   call compare("T spec", other)

   ! Pressure, with initial alpha and pressure
   other = critical_point(model, z0, zi, spec_CP%P, S=log(cp%P), max_iters=1000, a0=0.6_pr, P0=60.0_pr)
   call compare("P spec", other)

   ! Volume, with every initial guess
   other = critical_point(model, z0, zi, spec_CP%V, S=log(cp%Vx), max_iters=1000, &
      a0=0.6_pr, T0=400.0_pr, V0=0.2_pr, P0=60.0_pr)
   call compare("V spec", other)

   call test_critical_line()
   call test_critical_constants()

contains

   subroutine check_criticality(crit, name)
      !! The minimum eigenvalue of the stability matrix and the full set of
      !! criticality residuals vanish at the critical point.
      type(EquilibriumState), intent(in) :: crit
      character(len=*), intent(in) :: name
      real(pr) :: X(4)
      real(pr) :: F(4)
      real(pr) :: u(nc)
      real(pr) :: u_new(nc)
      real(pr) :: l1
      real(pr) :: P

      X = [log(crit%x(2)), log(crit%Vx), log(crit%T), log(crit%P)]
      u = [1.0_pr, 0.0_pr]
      l1 = lambda1(model, X, 0.0_pr, z0, zi, u, u_new=u_new, P=P)
      call assert(abs(l1) < 1e-6_pr, name // ": minimum eigenvalue of the stability matrix vanishes")
      call assert(abs(sqrt(sum(u_new**2)) - 1) < 1e-12_pr, name // ": normalized critical eigenvector")
      call assert(abs(P - crit%P) < 1e-8_pr*crit%P, name // ": lambda1 pressure equals the critical pressure")

      u = u_new
      F = F_critical(model, X, spec_CP%a, X(1), z0, zi, u)
      call assert(maxval(abs(F)) < 1e-5_pr, name // ": criticality residuals vanish")
   end subroutine check_criticality

   subroutine compare(name, crit)
      character(len=*), intent(in) :: name
      type(EquilibriumState), intent(in) :: crit
      call assert(crit%kind == "critical", name // ": converges")
      call assert(abs(crit%x(2) - a_spec) < 1e-6_pr, name // ": recovers the composition")
      call assert(abs(crit%T - cp%T) < 1e-5_pr*cp%T, name // ": same critical temperature")
      call assert(abs(crit%P - cp%P) < 1e-5_pr*cp%P, name // ": same critical pressure")
      call assert(abs(crit%Vx - cp%Vx) < 1e-5_pr*cp%Vx, name // ": same critical volume")
   end subroutine compare

   subroutine compare_with_envelope(crit)
      !! The PT envelope of the critical composition passes through the
      !! critical point.
      type(EquilibriumState), intent(in) :: crit
      type(EquilibriumState) :: dew
      type(PTEnvel2) :: env

      dew = saturation_temperature(model, crit%x, P=1.0_pr, kind="dew")
      env = pt_envelope_2ph(model, crit%x, dew)
      call assert(size(env%cps) >= 1, "PT envelope of the critical composition has a critical point")
      call assert(abs(env%cps(1)%T - crit%T) < 0.05_pr, "critical temperature matches the PT envelope")
      call assert(abs(env%cps(1)%P - crit%P) < 0.05_pr, "critical pressure matches the PT envelope")
   end subroutine compare_with_envelope

   subroutine test_critical_line()
      !! Type I critical line from methane to n-butane.
      type(CriticalLine) :: cl
      type(EquilibriumState) :: crit
      real(pr), parameter :: maxP = 100.0_pr
      integer :: n
      integer :: k

      ! With a pressure limit: stops before exceeding it
      cl = critical_line(model, a0=0.01_pr, z0=z0, zi=zi, ns0=spec_CP%a, S0=0.01_pr, dS0=0.01_pr, maxP=maxP)
      n = size(cl%P)
      call assert(n > 10, "critical line up to 100 bar has points")
      call assert(all(cl%P <= maxP), "critical line respects the maximum pressure")
      call assert(cl%P(n) > 0.9_pr*maxP, "critical line stops close to the maximum pressure")
      call assert(all(cl%z0 == z0) .and. all(cl%zi == zi), "critical line stores both fluids")
      call assert(size(cl%T) == n .and. size(cl%V) == n .and. size(cl%a) == n &
         .and. size(cl%ns) == n .and. size(cl%iters) == n, "critical line arrays have the same size")

      ! Every traced point is a critical point
      k = n/2
      crit = critical_point(model, z0, zi, spec_CP%a, S=cl%a(k), max_iters=1000, &
         a0=cl%a(k), T0=cl%T(k), V0=cl%V(k), P0=cl%P(k))
      call assert(abs(crit%T - cl%T(k)) < 1e-4_pr*cl%T(k), "critical line point: same T as critical_point")
      call assert(abs(crit%P - cl%P(k)) < 1e-4_pr*cl%P(k), "critical line point: same P as critical_point")

      ! Without a limit: a type I line ends at the n-butane critical point
      cl = critical_line(model, a0=0.01_pr, z0=z0, zi=zi, ns0=spec_CP%a, S0=0.01_pr, dS0=0.01_pr)
      n = size(cl%P)
      ! Where the tracing stops close to the pure end depends on rounding
      ! (alpha = 0.974 to 0.997 across builds), so check the approach only.
      call assert(cl%a(n) > 0.95_pr, "full critical line approaches pure n-butane")
      call assert(abs(cl%T(n) - 425.12_pr) < 5.0_pr, "full critical line approaches the n-butane Tc")
      call assert(abs(cl%P(n) - 37.96_pr) < 5.0_pr, "full critical line approaches the n-butane Pc")
      call assert(abs(cl%a(1) - 0.01_pr) < 1e-10_pr, "full critical line starts at the specified alpha")
      call assert(cl%T(1) > 190.564_pr .and. cl%T(1) < 200.0_pr, &
         "1% n-butane raises the critical temperature slightly above the methane Tc")
      call assert(maxval(cl%P) > 120.0_pr, "methane/n-butane critical line has a pressure maximum")
   end subroutine test_critical_line

   subroutine test_critical_constants()
      !! A cubic EoS built from (Tc, Pc, w) has its pure-component critical
      !! points exactly at (Tc, Pc), so these are a fixed point of
      !! `get_critical_constants`; the acentric factor from Psat(0.7 Tc) is
      !! close to the input one.
      type(CubicEoS) :: m3
      real(pr) :: tc(3)
      real(pr) :: pc(3)
      real(pr) :: w(3)

      tc = [190.564_pr, 425.12_pr, 304.21_pr]
      pc = [45.99_pr, 37.96_pr, 73.83_pr]
      w = [0.0115478_pr, 0.200164_pr, 0.223621_pr]
      m3 = PengRobinson76(tc, pc, w)
      ! The acentric factor is not used by the residual Helmholtz energy
      ! once the model is built, so it must be recomputed from Psat.
      m3%components%w = 0

      call get_critical_constants(m3)
      call assert(maxval(abs(m3%components%Tc - tc)/tc) < 1e-5_pr, "get_critical_constants recovers Tc")
      call assert(maxval(abs(m3%components%Pc - pc)/pc) < 1e-5_pr, "get_critical_constants recovers Pc")
      call assert(maxval(abs(m3%components%w - w)) < 5e-3_pr, "get_critical_constants recovers the acentric factor")
   end subroutine test_critical_constants
end program test_critical_point_specs
