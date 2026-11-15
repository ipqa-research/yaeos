module yaeos__equilibria_saturation_points
   !! # Saturation points
   !!
   !! # Description
   !! Bubble, dew, and liquid-liquid saturation points at specified
   !! temperature or pressure. The module first iterates a scalar Newton
   !! method on \(\sum_i z_i K_i - 1 = 0\) with successive-substitution
   !! updates of \(K_i\), and falls back to a full Newton solve in
   !! \(X = [\ln K_i, \ln T, \ln P]\) if that does not converge.
   !! A non-converged result is returned with `kind = "failed"`.
   !!
   !! # Examples
   !! ```fortran
   !! sat = saturation_pressure(model, z, T, kind="bubble")
   !! sat = saturation_temperature(model, z, P, kind="dew")
   !! ```
   use yaeos__constants, only: pr, R
   use yaeos__models_ar, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_auxiliar, only: k_wilson, P_wilson
   use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite

   implicit none

   real(pr) :: tol = 1e-6_pr
   !! Convergence tolerance for the saturation solvers
   integer :: max_iterations = 2000
   !! Maximum number of iterations for the full Newton solver
   integer :: iters_first_step = 30
   !! Maximum iterations for the initial scalar Newton step of saturation_pressure.
   integer :: iters_first_step_temperature = 100
   !! Maximum iterations for the initial scalar Newton step of saturation_temperature.
   real(pr) :: step_tol = 0.1_pr
   !! Currently unused (the scalar pressure step is limited to \(0.1P\)
   !! in code).

contains

   type(EquilibriumState) function saturation_pressure(model, n, t, kind, p0, y0, max_iters)
      !! # saturation_pressure
      !!
      !! Saturation pressure calculation function.
      !!
      !! ## Description
      !! Calculates the saturation pressure of a multicomponent mixture with
      !! a given molar composition `n`.
      !! It is possible to calculate:
      !!
      !! - Bubble point: `kind="bubble"`
      !! - Dew point: `kind="dew"`
      !! - Liquid-Liquid point: `kind="liquid-liquid"`
      !!
      !! It will first try to converge a solution using a \(1D\) Newton method to
      !! solve the equation
      !! \[
      !!    f(P) = \sum_i z_i K_i - 1 = 0
      !! \]
      !!
      !! updating \(K_i\) at each step as the ratio of fugacities of the phases.
      !! If the solution does not converge, it will use a full Newton method to
      !! solve the system in \(X = [\ln K_i, \ln T, \ln P]\) with \(\ln T\)
      !! specified. On failure the returned state has `kind = "failed"`.
      use yaeos__auxiliar, only: optval
      use yaeos__m_s_sp, only: solve_TP
      class(ArModel), target, intent(in) :: model !! EOS used to evaluate phase fugacities.
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar], defaults to the Wilson estimate
      real(pr), optional, intent(in) :: y0(:) !! Initial guess used as \(y\) to build \(K\) (`y/z` for bubble and liquid-liquid, `z/y` for dew); defaults to \(z K^{Wilson}\).
      integer, optional, intent(in) :: max_iters !! Currently ignored: the full Newton solve uses the module variable `max_iterations`.

      real(pr) :: P

      real(pr) :: k(size(n))
      real(pr) :: y(size(n)) !! Incipient/second phase composition.
      real(pr) :: z(size(n)) !! Feed or overall composition.
      real(pr) :: lnk(size(n))
      real(pr) :: lnfug_y(size(n))
      real(pr) :: dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n))
      real(pr) :: dlnphi_dp_z(size(n))
      real(pr) :: Vz
      real(pr) :: Vy

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f
      real(pr) :: step
      integer :: its !! Iteration count.
      integer :: iterations

      real(pr) :: X(size(n)+2)
      real(pr) :: S
      integer :: ns
      integer :: nc !! Number of components.

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (present (p0)) then
         p = p0
      else
         P = p_wilson(model, z, T)
      end if

      if (present(y0)) then
         y = y0
      else
         y = z * k_wilson(model, T, P)
      end if
      iterations = optval(max_iters, max_iterations)

      select case(kind)
       case("bubble")
         k = y/z
         incipient = "vapor"
         main = "liquid"
       case("dew")
         k = z/y
         incipient = "liquid"
         main = "vapor"
       case("liquid-liquid")
         k = y/z
         incipient = "liquid"
         main = "liquid"
      end select

      where (z == 0)
         k = 0
      end where
      ! ========================================================================

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      do its=1, iters_first_step
         y = k*z
         call model%lnphi_pt(y, P, T, vy, incipient, lnPhi=lnfug_y, dlnphidp=dlnphi_dp_y)
         call model%lnphi_pt(z, P, T, vz, main, lnPhi=lnfug_z, dlnphidp=dlnphi_dp_z)

         k = exp(lnfug_z - lnfug_y)

         if (all(k < 1e-9_pr) .or. all(abs(k-1) < tol)) exit

         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))

         do while (P - step < 0 .or. abs(step) > 0.1*P)
            step = step/2
         end do

         p = p - step
         if (abs(step) < tol .and. abs(f) < tol) exit
      end do
      ! ========================================================================
      if (its >= iters_first_step) then
         nc = size(n)
         X(:nc) = log(y/z)
         X(nc+1) = log(T)
         X(nc+2) = log(P)
         ns = nc+1
         S = X(ns)

         call solve_TP(model, kind, z, X, ns, S, tol, max_iterations, its)

         P = exp(X(nc+2))
         y = z * exp(X(:nc))
         call model%volume(n=z, P=P, T=T, V=Vz, root_type=main)
         call model%volume(n=y, P=P, T=T, V=Vy, root_type=incipient)
      end if

      y = y/sum(y)

      select case(kind)
       case("bubble")
         saturation_pressure = EquilibriumState(kind="bubble", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
       case("dew")
         saturation_pressure = EquilibriumState(kind="dew", &
            iters=its, x=y, y=z, vy=vz, vx=vy, t=t, p=p, beta=1._pr&
            )
       case("liquid-liquid")
         saturation_pressure = EquilibriumState(kind="liquid-liquid", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
      end select

      if (its >= max_iterations .or. any(ieee_is_nan(y)) .or. ieee_is_nan(P) .or. abs(Vz - Vy) < 1e-8) then
         saturation_pressure%kind = "failed"
      end if
   end function saturation_pressure

   type(EquilibriumState) function saturation_temperature(model, n, p, kind, t0, y0, max_iters)
      !! # saturation_temperature
      !!
      !! # Description
      !! Saturation temperature calculation function.
      !!
      !! Calculates the saturation temperature of a multicomponent mixture with
      !! a given molar composition `n`.
      !! It is possible to calculate:
      !!
      !! - Bubble point: `kind="bubble"`
      !! - Dew point: `kind="dew"`
      !! - Liquid-Liquid point: `kind="liquid-liquid"`
      !!
      !! It will first try to converge a solution using a \(1D\) Newton method to
      !! solve the equation
      !! \[
      !!    f(T) = \sum_i z_i K_i - 1 = 0
      !! \]
      !!
      !! updating \(K_i\) at each step as the ratio of fugacities of the phases.
      !! If the solution does not converge, it will use a full Newton method to
      !! solve the system in \(X = [\ln K_i, \ln T, \ln P]\) with \(\ln P\)
      !! specified. On failure the returned state has `kind = "failed"`.
      use yaeos__auxiliar, only: optval
      use yaeos__m_s_sp, only: solve_TP
      class(ArModel), target, intent(in) :: model !! EOS used to evaluate phase fugacities.
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Pressure [bar]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: t0 !! Initial temperature [K], defaults to 250 K
      real(pr), optional, intent(in) :: y0(:) !! Initial guess used as \(y\) to build \(K\) (`y/z` for bubble and liquid-liquid, `z/y` for dew); defaults to \(z K^{Wilson}\).
      integer, optional, intent(in) :: max_iters !! Currently ignored: the full Newton solve uses the module variable `max_iterations`.

      real(pr) :: t
      real(pr) :: vy
      real(pr) :: vz

      real(pr) :: k(size(n))
      real(pr) :: y(size(n)) !! Incipient/second phase composition.
      real(pr) :: z(size(n)) !! Feed or overall composition.
      real(pr) :: lnk(size(n))
      real(pr) :: lnfug_y(size(n))
      real(pr) :: dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n))
      real(pr) :: dlnphi_dt_z(size(n))

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f
      real(pr) :: step
      integer :: its !! Iteration count.
      integer :: iterations

      logical :: is_incipient(size(n)) !! `.false.` for components absent from the initial `y`.

      real(pr) :: T1
      real(pr) :: T2
      real(pr) :: f1 !! Unused.
      real(pr) :: f2 !! Unused.

      integer :: nc !! Number of components.
      integer :: ns
      real(pr) :: X(size(n)+2)
      real(pr) :: S

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      is_incipient = .true.
      z = n/sum(n)
      if (present (t0)) then
         t = t0
      else
         t = 250._pr
      end if

      if (present(y0)) then
         y = y0
      else
         y = z * k_wilson(model, T, P)
         y = y/sum(y)
      end if
      iterations = optval(max_iters, max_iterations)

      select case(kind)
       case("bubble")
         k = y/z
         incipient = "vapor"
         main = "liquid"
       case("dew")
         k = z/y
         incipient = "liquid"
         main = "vapor"
       case("liquid-liquid")
         k = y/z
         incipient = "liquid"
         main = "liquid"
      end select

      where (z == 0)
         k = 0
      end where

      where (y == 0)
         is_incipient = .false.
      end where

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      f = 10
      its = 0
      do while(abs(f) > tol .and. its < iters_first_step_temperature)
         its = its + 1
         y = k*z
         call model%lnphi_pt(y, P, T, vy, incipient, lnPhi=lnfug_y, dlnphidt=dlnphi_dt_y)
         call model%lnphi_pt(z, P, T, vz, main, lnPhi=lnfug_z, dlnphidt=dlnphi_dt_z)

         k = exp(lnfug_z - lnfug_y)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))

         if (T - step < 50) then
            T = 100 + abs(step)
         end if

         do while(abs(step) > 50)
            step = step/2
         end do

         T = T - step
         if (abs(step) < tol .and. abs(f) < tol) exit
      end do

      if (its >= iters_first_step_temperature) then
         nc = size(n)
         X(:nc) = log(k)
         X(nc+1) = log(T)
         X(nc+2) = log(P)
         ns = nc+2
         S = X(ns)

         call solve_TP(model, kind, z, X, ns, S, tol, max_iterations, its)

         T = exp(X(nc+1))
         y = z * exp(X(:nc))
         call model%volume(n=z, P=P, T=T, V=Vz, root_type=main)
         call model%volume(n=y, P=P, T=T, V=Vy, root_type=incipient)
      end if

      select case(kind)
       case("bubble")
         saturation_temperature = EquilibriumState(kind="bubble", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
       case("dew")
         saturation_temperature = EquilibriumState(kind="dew", &
            iters=its, x=y, y=z, vy=vz, vx=vy, t=t, p=p, beta=1._pr&
            )
       case("liquid-liquid")
         saturation_temperature = EquilibriumState(kind="liquid-liquid", &
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
            )
      end select

      if (its >= max_iterations .or. any(ieee_is_nan(y)) .or. ieee_is_nan(T) .or. abs(Vz - Vy) < 1e-8) then
         saturation_temperature%kind = "failed"
      end if
   end function saturation_temperature

end module yaeos__equilibria_saturation_points
