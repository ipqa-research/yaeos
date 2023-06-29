module yaeos_thermo_properties
   !-> Thermodynamic routines
   ! This module includes all the procedures to calculate bulk properties like
   ! Pressure, fugacity, volume, etc.
   !
   ! This procedures are based on residual Helmholtz free energy models. Using
   ! The approach presented by Michelsen and Møllerup.
   use constants, only: pr, R
   use yaeos_models, only: residual_helmholtz
   use hyperdual_mod
   use yaeos_interfaces, only: volume_initalizer
   implicit none

   private

   public :: pressure, get_volume, ln_phi, vinit

   procedure(volume_initalizer), pointer :: vinit

contains

   ! =============================================================================
   !  Bulk Properties
   ! -----------------------------------------------------------------------------
   subroutine pressure(z, v, t, p, dp, dp2)
      !-> Calculate pressure using the residual Helmholtz energy, as defined
      !   by Michelsen Møllerup.
      real(pr), intent(in)            :: z(:) !- Molar compositions
      real(pr), intent(in)            :: v    !- Volume
      real(pr), intent(in)            :: t    !- Temperature

      real(pr), intent(out)           :: p    !- Pressure
      real(pr), optional, intent(out) :: dp(size(z) + 2) !- Pressure derivatives
      real(pr), optional, intent(out) :: dp2(size(z) + 2, size(z) + 2)

      ! Ar values and derivatives
      real(pr) :: ar, dar(size(z) + 2), dar2(size(z) + 2, size(z) + 2)

      ! Keep the RT product
      real(pr) :: RT

      integer :: n

      n = size(z)
      RT = R*t

      call residual_helmholtz(z, v, t, ar, dar, dar2)

      p = -RT*dar(n + 1) + sum(z)*RT/v
      if (present(dp)) then
         dp(:n) = -RT*dar2(:n, n + 1) + Rt/v
         dp(n + 1) = -RT*dar2(n + 1, n + 1) - sum(z)*RT/v**2
      end if
      ! dp(n+2) = -RT * dar2(n+1, n+2) + p/t
   end subroutine

   subroutine ln_phi(z, v, t, lnphi)
      !-> Fugacity coefficent as a function of volume and temperature.
      !   The fugacity coefficent is the first derivative of the Residual
      !   Helmholtz energy, substracted with the compressibility factor:
      !
      !   \[ln \phi_i = \frac{dAr}{dn_i}(z, v, T) - \ln Z\]
      !
      real(pr), intent(in)  :: z(:)           !- Composition vector
      real(pr), intent(in)  :: v              !- Volume
      real(pr), intent(in)  :: t              !- Temperature
      real(pr), intent(out) :: lnphi(size(z)) !- Fugacity coefficent

      real(pr) :: ar, dar(size(z) + 2), dar2(size(z) + 2, size(z) + 2)
      real(pr) :: p, compressibility

      call pressure(z, v, t, p, dar, dar2)
      call residual_helmholtz(z, v, t, ar, dar, dar2)

      compressibility = p*v/(R*t)
      lnphi = dar(:size(z)) - log(compressibility)
   end subroutine

   recursive subroutine get_volume(z, p, t, v, root, v0)
      !-> Volume solver routine.
      !   This volume solving routine uses a simple Newton method to solve the
      !   equation:
      !
      !   \[ 0 = P(z, v, T) - P_{spec} \]
      !
      !   I'ts expected from the model to include an initialization function
      !   for the first volume ( $v_0$ ), besides that, an initialization
      !   value can be provided as an optional argument.
      real(pr), intent(in)           :: z(:) !- Composition Vector
      real(pr), intent(in)           :: p    !- Pressure
      real(pr), intent(in)           :: t    !- Temperature
      character(len=*), intent(in)   :: root !- root: `[vapor | liquid | stable]`
      real(pr), intent(in), optional :: v0   !- Initial volume

      real(pr), intent(out)          :: v    !- Volume

      real(pr) :: p_in, dp(size(z) + 2), dp2(size(z) + 2, size(z) + 2), delta
      real(pr) :: v_liq, v_vap

      integer :: its, n
      n = size(z)

      if (present(v0)) then
         v = v0
      else
         select case (root)
         case ("vapor")
            v = R*t/p
         case ("liquid")
            v = 2.0_pr*vinit(z, p, t)
         case ("stable")
            call get_volume(z, v, t, v_liq, "liquid")
            call get_volume(z, v, t, v_vap, "vapor")
         end select
      end if

      delta = 1
      p_in = p*20
      its = 0

      do while ((abs(delta) > 1e-10 .or. abs(p - p_in) > 1e-10) &
                .and. its < 100)
         its = its + 1
         call pressure(z, v, t, p_in, dp, dp2)
         delta = ((p - p_in)/dp(n + 1)) ! / p
         v = v + delta
      end do
   end subroutine
   ! ===========================================================================
end module
