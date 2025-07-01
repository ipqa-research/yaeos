module yaeos__models_ge_nrtlhv
   !! NRTL-HV model for excess Gibbs energy
   use yaeos__constants, only: pr, R
   use yaeos__autodiff
   use yaeos__models_ge, only: GeModel

   implicit none

   type, extends(GeModel) :: NRTLHV
      !! # NRTLHV
      !! Huron-Vidal modification of the NRTL model for excess Gibbs energy.
      !!
      !! # Description
      !! Huron-Vidal modified the NRTL model, including the covolume parameter
      !! into the model. This helps when running the model coupled with a cubic
      !! equation of state. Because it can be reduced to the classical
      !! quadratic mixing rules when some parameters are defined accordingly
      !! # Examples
      !!
      !! # References
      real(pr), allocatable :: b(:) !! Covolume parameter
      real(pr), allocatable :: alpha(:, :) !! \( \alpha \) matrix
      real(pr), allocatable :: gij(:, :) !! \( g_{ij} \) matrix
   contains
      procedure :: excess_gibbs => excess_gibbs
   end type NRTLHV

contains

   subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Calculate Excess Gibbs and its derivatives.
      use yaeos__models_ge_base, only: nrtl_hv_ge, nrtl_hv_tdep
      class(NRTLHV), intent(in) :: self !! Model
      real(pr), intent(in) ::n(:) !! Moles vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs free energy
      real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), optional, intent(out) :: Gen(size(n)) !! \(\frac{dG^E}{dn}\)
      real(pr), optional, intent(out) :: GeTn(size(n))
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))
      real(pr) :: tau(size(n), size(n)), dtaudt(size(n), size(n)), dtaudt2(size(n), size(n))

      call nrtl_hv_tdep(T, self%gij, tau, dtaudt, dtaudt2)
      call nrtl_hv_ge(n=n, T=T,&
         b=self%b, alpha=self%alpha, &
         tau=tau, dtaudt=dtaudt, dtaudt2=dtaudt2, &
         Ge=Ge, Gen=Gen, GeT=GeT, GeT2=GeT2, GeTn=GeTn, Gen2=Gen2)
   end subroutine excess_gibbs


   ! subroutine C_from_cubic(self, n, T, C)
   !    class(NRTLHV), intent(in) :: self
   !    type(hyperdual), intent(in) :: n(:)
   !    type(hyperdual), intent(in) :: T
   !    type(hyperdual), intent(out) :: C

   !    real(pr) :: a(size(n)), dadt(size(n)), dadt2(size(n))
   !    real(pr) :: aij(size(n), size(n))
   !    real(pr) :: daijdt(size(n), size(n))
   !    real(pr) :: daijdt2(size(n), size(n))
   !    real(pr) :: Tr(size(n))
   !    real(pr) :: d1, dd1i(size(n)), dd1ij(size(n), size(n))
   !    real(pr) :: lambda, lambdadn(size(n)), dlambdadn2(size(n), size(n))
   !    integer :: i, j, nc

   !    type(hyperdual) :: g_ii(size(n)), g_ji(size(n), size(n))
   !    type(hyperdual) :: a_hd(size(n)), b_hd(size(n))
   !    nc = size(n)

   !    Tr = T%f0/self%components%Tc

   !    !
   ! ========================================================================
   !    ! Attractive parameter and derivatives
   !    !
   ! ------------------------------------------------------------------------
   !    call self%cubic%mixrule%D1mix(n%f0, self%cubic%del1, d1, dd1i, dd1ij)
   !    call lamdba_hv(d1, dd1i, dd1ij, lambda, lambdadn, dlambdadn2)
   !
   !    call    self%cubic%alpha%alpha(Tr, a, dadt, dadt2)
   !    a =     self%cubic%ac * a
   !    dadt =  self%cubic%ac * dadt / self%cubic%components%Tc
   !    dadt2 = self%cubic%ac * dadt2 / self%cubic%components%Tc**2

   !    g_ii = -a_hd/b_hd * lambda
   !    do i=1,nc
   !       do j=1,nc
   !          if (.not. self%use_cubic(i, j)) cycle
   !          g_ji(j, i) = -2*lambda * sqrt(a(i)*a(j)) * 1/(b(i)+b(j)) *
   ! (1-kij(i, j))
   !       end do
   !    end do

   !    do i=1,nc
   !       do j=1,nc
   !          C(j, i) = (g_ji(j, i) - g_ji(i,i))
   !       end do
   !    end do

   ! end subroutine C_from_cubic
end module yaeos__models_ge_nrtlhv
