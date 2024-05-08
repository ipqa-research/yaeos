module yaeos_equilibria_saturation_points
   use yaeos_constants, only: pr
   use yaeos_models, only: ArModel
   use yaeos_thermoprops, only: fugacity_vt, fugacity_tp
   use yaeos_equilibria_equilibria_state, only: EquilibriaState
   use yaeos__phase_equilibria_auxiliar, only: k_wilson

   real(pr) :: tol = 1e-9_pr
   integer :: max_iterations = 100
   real(pr) :: step_tol = 0.1_pr

contains

   type(EquilibriaState) function bubble_pressure(model, n, t, p0, y0, max_inner_its)
      !! Bubble pressure calculation function.
      !!
      !! Calculates the bubble temperature of a multicomponent mixture.
      use stdlib_optval, only: optval
      use yaeos_thermoprops, only: pressure
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: n(:) !! Composition vector [moles / molar fraction]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      real(pr), optional, intent(in) :: y0(:) !! Initial composition
      integer, optional, intent(in) :: max_inner_its(:) !! Inner iterations

      real(pr) :: p, vy, vz

      real(pr) :: k(size(n)), y(size(n)), z(size(n)), lnk(size(n))
      real(pr) :: lnfug_y(size(n)), dlnphi_dp_y(size(n))
      real(pr) :: lnfug_z(size(n)), dlnphi_dp_z(size(n))

      real(pr) :: f, step
      integer :: its, inner_its

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
      inner_its = optval(inner_its, 50)
      ! ========================================================================

      ! ========================================================================
      !  Solve point
      ! ------------------------------------------------------------------------
      do its=1,max_iterations
         call fugacity_tp(model, y, T, P, vy, "vapor", lnphip=lnfug_y, dlnphidp=dlnphi_dp_y)
         call fugacity_tp(model, z, T, P, vz, "liquid", lnphip=lnfug_z, dlnphidp=dlnphi_dp_z)

         do while (any(isnan(lnfug_y)))
            inner_its = inner_its + 1
            p = p/2.0_pr
            call fugacity_tp(model, y, T, P, vy, "vapor", lnphip=lnfug_y, dlnphidp=dlnphi_dp_y)
            call fugacity_tp(model, z, T, P, vz, "liquid", lnphip=lnfug_z, dlnphidp=dlnphi_dp_z)
         end do

         lnk = lnfug_z - lnfug_y
         k = exp(lnk)

         f = sum(z*k) - 1
         step = f/sum(z * k * (dlnphi_dp_z - dlnphi_dp_y))
         p = p - step
         y = z*k

         if (abs(step) < tol) exit
      end do
      bubble_pressure = EquilibriaState(&
         iters=its, y=y, x=z, vx=vz, vy=vy, t=t, p=p, beta=0._pr&
         )
      ! ========================================================================
   end function bubble_pressure
end module yaeos_equilibria_saturation_points
