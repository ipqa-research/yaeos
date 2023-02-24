module thermo_properties
   !! Thermodynamic routines
   use constants
   use models, only: ArModel, residual_helmholtz
   implicit none 

   interface pressure
      module procedure :: ArModel_pressure
   end interface pressure

contains

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

      p = -Rt*dar(n+1) + sum(z)*Rt / v
      dp(:n) = -Rt*dar2(:n, n+1) + Rt / v
      dp(n+1) = -Rt*dar2(n+1, n+1) - sum(z)*Rt/v**2
      dp(n+2) = -Rt * dar2(n+1, n+2) + p/t

   end subroutine

   ! subroutine get_volume(model, z, p, t, v)
   ! end subroutine

   subroutine lnphi()
   end subroutine lnphi

end module thermo_properties