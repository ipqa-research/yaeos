module yaeos__equilibria_saturation_points
   use yaeos__constants, only: pr, R
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_auxiliar, only: k_wilson, P_wilson
   use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite

   implicit none

   real(pr) :: tol = 1e-6_pr
   integer :: max_iterations = 2000
   integer :: iters_first_step = 15
   real(pr) :: step_tol = 0.1_pr

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
      !! solve the system of equations using the variables \(K_i\) and \(\ln P\).
      use yaeos__auxiliar, only: optval
      use yaeos__m_s_sp, only: solve_TP
      class(ArModel), target, intent(in) :: model
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_iters !! Maximum number of iterations

      real(pr) :: P

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))
      real(pr) :: Vz, Vy

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f, step
      integer :: its, iterations

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
         block
            real(pr) :: X(size(n)+2), S
            integer :: ns, nc
            nc = size(n)
            X(:nc) = log(y/z)
            X(nc+1) = log(T)
            X(nc+2) = log(P)
            ns = nc+1
            S = X(ns)

            call solve_TP(model, kind, z, X, ns, S, tol, max_iterations, its)

            P = exp(X(nc+2))
            y = z * exp(X(:nc))
            call model%volume(n=n, P=P, T=T, V=Vz, root_type=main)
            call model%volume(n=y, P=P, T=T, V=Vy, root_type=incipient)
         end block
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

      if (its >= max_iterations .or. any(isnan([P, y])) .or. abs(Vz - Vy) < 1e-8) then
         saturation_pressure%kind = "failed"
      end if
   end function saturation_pressure

   type(EquilibriumState) function saturation_temperature(model, n, p, kind, t0, y0, max_iters)
      !! Saturation temperature calculation function.
      !!
      !! Calculates the saturation pressure of a multicomponent mixture with
      !! a given molar composition `n`.
      !! It is possible to calculate:
      !!
      !! - Bubble point: `kind="bubble"`
      !! - Dew point: `kind="dew"`
      !! - Liquid-Liquid point: `kind="liquid-liquid"`
      !! It will first try to converge a solution using a \(1D\) Newton method to
      !! solve the equation
      !! \[
      !!    f(P) = \sum_i z_i K_i - 1 = 0
      !! \]
      !!
      !! updating \(K_i\) at each step as the ratio of fugacities of the phases.
      !! If the solution does not converge, it will use a full Newton method to
      !! solve the system of equations using the variables \(K_i\) and \(\ln T\).
      use yaeos__auxiliar, only: optval
      use yaeos__m_s_sp, only: solve_TP
      class(ArModel), target, intent(in) :: model
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: p !! Pressure [bar]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: t0 !! Initial temperature [K]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_iters !! Maximum number of iterations

      real(pr) :: t, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dt_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dt_z(size(n))

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f, step
      integer :: its, iterations

      logical :: is_incipient(size(n))

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
      do its=1, 5
         y = k*z
         where (.not. is_incipient)
            y = 0
         endwhere

         call model%lnphi_pt(y, P, T, vy, incipient, lnPhi=lnfug_y, dlnphidt=dlnphi_dt_y)
         call model%lnphi_pt(z, P, T, vz, main, lnPhi=lnfug_z, dlnphidt=dlnphi_dt_z)

         k = exp(lnfug_z - lnfug_y)
         f = sum(z*k) - 1
         step = f/sum(T * z * k * (dlnphi_dt_z - dlnphi_dt_y))

         if (.not. ieee_is_finite(step) .or. ieee_is_nan(step)) exit

         do while (T - step < 0)
            if (isnan(step)) step = 10
            step = step/2
         end do

         t = t - step

         if (abs(step) < tol .and. abs(f) < tol) exit
      end do
      ! ========================================================================
      its = iters_first_step
      if (its >= iters_first_step) then
         block
            real(pr) :: X(size(n)+2), S
            integer :: ns, nc
            nc = size(n)
            X(:nc) = log(y/z)
            X(nc+1) = log(T)
            X(nc+2) = log(P)
            ns = nc+2
            S = X(ns)

            call solve_TP(model, kind, z, X, ns, S, tol, max_iterations, its)

            T = exp(X(nc+1))
            y = z * exp(X(:nc))
            call model%volume(n=n, P=P, T=T, V=Vz, root_type=main)
            call model%volume(n=y, P=P, T=T, V=Vy, root_type=incipient)
         end block
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

      if (its >= max_iterations .or. any(isnan([P, y])) .or. abs(Vz - Vy) < 1e-8) then
         saturation_temperature%kind = "failed"
      end if
   end function saturation_temperature

end module yaeos__equilibria_saturation_points
