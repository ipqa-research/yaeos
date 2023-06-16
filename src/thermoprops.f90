module thermo_properties
   !! Thermodynamic routines
   use constants
   use models, only: ArModel, residual_helmholtz
   use hyperdual_mod
   implicit none

   private

   public :: pressure, volume, ln_phi

   interface volume
      module procedure :: ArModel_get_volume !! Volume solving routine
   end interface

   interface ln_phi
      !! Fugacity coefficent
      module procedure :: ArModel_ln_phi
   end interface

   interface pressure
         module procedure :: ArModel_pressure
   end interface
contains

   ! =============================================================================
   !  Bulk Properties
   ! -----------------------------------------------------------------------------
   subroutine ArModel_pressure(model, z, v, t, p, dp, dp2)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(model%size)
      real(pr), intent(in) :: v, t

      real(pr), intent(out) :: p
      real(pr), intent(out) :: dp(model%size + 2)
      real(pr), intent(out) :: dp2(model%size + 2, model%size + 2)

      real(pr) :: ar, dar(model%size + 2), dar2(model%size + 2, model%size + 2)
      real(pr) :: RT

      integer :: n
      
      n = model%size
      RT = R*t

      call residual_helmholtz(model, z, v, t, ar, dar, dar2)

      p = -RT*dar(n+1) + sum(z)*RT / v
      dp(:n)  = -RT * dar2(:n, n+1) + Rt / v
      dp(n+1) = -RT * dar2(n+1, n+1) - sum(z)*RT/v**2
      ! dp(n+2) = -Rt * dar2(n+1, n+2) + p/t
   end subroutine

   subroutine ArModel_ln_phi(model, z, v, t, lnphi)
      !! Fugacity coefficent as a function of volume and temperature.
      !! The fugacity coefficent is the first derivative of the Residual
      !! Helmholtz energy, substracted with the compressibility factor:
      !!
      !! \[ln \phi_i = \frac{dAr}{dn_i}(z, v, T) - \ln Z\]
      !!
      class(ArModel), intent(in) :: model !! Residual Helmholtz Mode Derived type
      real(pr), intent(in) :: z(model%size) !! Composition vector
      real(pr), intent(in) :: v !! Volume
      real(pr), intent(in) :: t !! Temperature
      real(pr), intent(out) :: lnphi(model%size) !! Fugacity coefficent

      real(pr) :: ar, dar(model%size + 2), dar2(model%size + 2, model%size + 2)
      real(pr) :: p, compressibility

      call pressure(model, z, v, t, p, dar, dar2)
      call residual_helmholtz(model, z, v, t, ar, dar, dar2)

      compressibility = p * v / (R * t)
      lnphi = dar(:model%size) - log(compressibility)
   end subroutine

   subroutine ArModel_get_volume(model, z, p, t, v, root, v0)
      !! Volume solver routine.
      !! This volume solving routine uses a simple Newton method to solve the
      !! equation:
      !!
      !! \[ 0 = P(z, v, T) - P_{spec} \]
      !!
      !! I'ts expected from the model to include an initialization function 
      !! for the first volume ( $v_0$ ), besides that, an initialization 
      !! value can be provided as an optional argument.
      !!
      ! TODO: The optional should be a function with no parameters instead
      class(ArModel), intent(in) :: model !! Residual Helmholtz Model
      real(pr), intent(in) :: z(model%size) !! Composition Vector
      real(pr), intent(in) :: p !! Pressure
      real(pr), intent(in) :: t !! Temperature
      character(len=*), intent(in) :: root  !! Desired root `[vapor | liquid | stable]`
      real(pr), intent(in), optional :: v0 !! Initial volume
      real(pr), intent(out) :: v !! Volume

      real(pr) :: p_in, dp(model%size + 2), dp2(model%size + 2, model%size + 2), delta

      integer :: it, n
      n = model%size

      if (present(v0)) then
         v = v0
      else
         select case (root)
         case ("vapor")
            v = R * t / p
         case ("liquid")
            ! TODO: Generalize the initialization with a function
            v = 0.03
         case ("stable")
         end select
      end if

      delta = 1

      do while (abs(delta) > 1e-10)
         call pressure(model, z, v, t, p_in, dp, dp2)
         delta = (p - p_in)/dp(n + 1)
         v = v + delta
      end do
   end subroutine
   ! =============================================================================
end module thermo_properties