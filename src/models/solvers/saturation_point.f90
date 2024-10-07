module yaeos__m_s_sp
   !! Module to calculate saturation points
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel, size
   implicit none

contains

   subroutine saturation_F(model, X, ns, S, F)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(out) :: F(:)

      ! Variables
      real(pr) :: T, V_main, V_inc
      real(pr) :: z(size(model))


      ! Incipient phase variables
      real(pr) :: y(size(z))
      real(pr) :: lnphi_inc(size(model))

      real(pr) :: lnphi_main(size(model))
      real(pr) :: P_main, P_inc

      integer :: nc

      nc = size(z)

      y = z * exp(lnK)
      
      call model%lnphi_vt(z, V_main, T, P_main, lnPhi=lnphi_main)
      call model%lnphi_vt(y, V_inc, T, P_inc, lnPhi=lnphi_inc)

      F(1:nc) = lnK - lnphi_main + lnphi_inc
      F(nc + 1) = sum(y - z)
      F(nc + 2) = P_main - P_inc
   end subroutine
end module