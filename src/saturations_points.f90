module yaeos__equilibria_saturation_points
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_auxiliar, only: k_wilson

   real(pr) :: tol = 1e-9_pr
   integer :: max_iterations = 1000
   real(pr) :: step_tol = 0.1_pr

contains

   type(EquilibriumState) function saturation_pressure(model, n, t, kind, p0, y0, max_iters)
      !! Saturation pressure calculation function.
      !!
      !! Calculates the saturation pressure of a multicomponent mixture with
      !! a given molar composition `n`.
      !! It is possible to calculate:
      !!
      !! - Bubble point: `kind="bubble"`
      !! - Dew point: `kind="dew"`
      !! - Liquid-Liquid point: `kind="liquid-liquid"`
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      character(len=*), intent(in) :: kind !! [bubble|dew|liquid-liquid]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_iters !! Maximum number of iterations

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))

      character(len=50) :: incipient
      character(len=50) :: main

      real(pr) :: f, step
      integer :: its, iterations, i

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (present (p0)) then
         p = p0
      else
         call model%pressure(z, T, 10._pr, P=P)
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
      do its=1, iterations
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
      use stdlib_optval, only: optval
      class(ArModel), intent(in) :: model
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
         t = 150._pr
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

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      do its=1, iterations
         y = k*z
         where (.not. is_incipient)
            y = 0
         endwhere

         call model%lnphi_pt(y, P, T, vy, incipient, lnPhi=lnfug_y, dlnphidt=dlnphi_dt_y)
         call model%lnphi_pt(z, P, T, vz, main, lnPhi=lnfug_z, dlnphidt=dlnphi_dt_z)

         k = exp(lnfug_z - lnfug_y)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dt_z - dlnphi_dt_y))

         do while (abs(step) > 0.25*T .or. T - step < 0)
            step = step/2
         end do

         t = t - step

         if (abs(step) < tol .and. abs(f) < tol) exit
      end do
      ! ========================================================================
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
   end function saturation_temperature
end module yaeos__equilibria_saturation_points
