module yaeos_equilibria_saturation_points
   use yaeos_constants, only: pr
   use yaeos_models, only: ArModel
   use yaeos_thermoprops, only: fugacity_vt, fugacity_tp
   use yaeos_equilibria_equilibria_state, only: EquilibriaState
   use yaeos__phase_equilibria_auxiliar, only: k_wilson

   real(pr) :: tol = 1e-9_pr
   integer :: max_iterations = 1000
   real(pr) :: step_tol = 0.1_pr

contains

   type(EquilibriaState) function saturation_pressure(model, n, t, kind, p0, y0, max_iters)
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
      use yaeos_thermoprops, only: pressure
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
      integer :: its, iterations

      ! =======================================================================
      ! Handle arguments
      ! -----------------------------------------------------------------------
      z = n/sum(n)
      if (present (p0)) then
         p = p0
      else
         call pressure(model, z, T, 10._pr, P=P)
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
         incipient = "liquid"
         main = "liquid"
      end select
      ! ========================================================================

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      do its=1, iterations
         y = k*z
         call fugacity_tp(model, y, T, P, vy, incipient, lnphip=lnfug_y, dlnphidp=dlnphi_dp_y)
         call fugacity_tp(model, z, T, P, vz, main, lnphip=lnfug_z, dlnphidp=dlnphi_dp_z)

         k = exp(lnfug_z - lnfug_y)
         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))

         do while (abs(step) > 0.1*P)
            step = step/2
         end do

         p = p - step
         if (abs(step) < tol .and. abs(f) < tol) exit
      end do
      ! ========================================================================
      select case(kind)
       case("bubble")
         saturation_pressure = EquilibriaState(&
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
         )
       case("dew")
         saturation_pressure = EquilibriaState(&
            iters=its, x=y, y=z, vy=vz, vx=vy, t=t, p=p, beta=1._pr&
       )
       case("liquid-liquid")
         saturation_pressure = EquilibriaState(&
            iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
      )
      end select
   end function saturation_pressure
end module yaeos_equilibria_saturation_points
