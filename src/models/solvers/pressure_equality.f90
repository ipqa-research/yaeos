module yaeos__solvers_pressure_equality
   !! Solve the pressure equality of a
   use yaeos__constants, only: pr, R
   use yaeos__models_ar, only: ArModel

   implicit none

contains

   subroutine pressure_equality_V_beta_xy(model, T, V, beta, x, y, vx, vy, P)
      !! Solve pressure equality between two phases at a given temperature,
      !! total volume, vapor molar fractions and compositions.
      use iso_fortran_env, only: error_unit

      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: V !! Total volume [L/mol]
      real(pr), intent(in) :: beta !! Molar fraction of light-phase
      real(pr), intent(in) :: x(:) !! Molar fractions of heavy-phase
      real(pr), intent(in) :: y(:) !! Molar fractions of light-phase
      real(pr), intent(in out) :: Vx !! Heavy-phase molar volume [L/mol]
      real(pr), intent(in out) :: Vy !! Light-Phase molar volume [L/mol]
      real(pr), intent(out) :: P !! Pressure [bar]

      real(pr) :: Bx !! Liquid phase covolume
      real(pr) :: dVydVx !! Derivative of Vy wrt Vx

      ! Pressure equality newton functions
      real(pr) :: h !! Pressure equality
      real(pr) :: dh  !! dh/
      real(pr) :: stepv

      real(pr) :: dPxdV, dPydV
      real(pr) :: Px, Py

      integer :: its

      dVydVx = -(1 - beta)/beta
      Bx = model%get_v0(x, 0.1_pr, T)

      ! First evaluation will be with Vx = 1.5*Bx
      if (Vx < Bx) Vx = 1.625_pr*Bx

      call model%pressure(x, Vx, T, Px, dpdv=dPxdV)

      do while (Px < 0 .or. dPxdV >= 0)
         Vx = Vx - 0.2*(Vx - Bx)
         call model%pressure(x, Vx, T, Px, dpdv=dPxdV)
      end do

      Vy = (V - (1 - beta)*Vx)/beta

      h = 1.0
      its = 0
      do while (abs(h) > 1.d-4)
         ! Newton for solving P equality, with Vx as independent variable
         its = its + 1

         call model%pressure(x, Vx, T, Px, dpdv=dPxdV)
         call model%pressure(y, Vy, T, Py, dpdv=dPydV)

         h = Py - Px
         dh = -dPydV * dVydVx - dPxdV
         stepv = -h/dh

         if (its >= 10) stepv = stepv/2

         Vx = Vx + stepv

         do while (Vx < 1.001*Bx)
            stepv = stepv/2
            Vx = Vx - stepv
         end do

         Vy = (v - (1 - beta)*Vx)/beta

         if (its >= 100) then
            write (error_unit, *) "WARN(FLASH_VT): volume convergence problems", Px, Py
            P = -1.0
            return
         end if
      end do

      call model%pressure(x, Vx, T, Px)
      call model%pressure(y, Vy, T, Py)
      P = (Px + Py) * 0.5_pr
   end subroutine pressure_equality_V_beta_xy
end module yaeos__solvers_pressure_equality
