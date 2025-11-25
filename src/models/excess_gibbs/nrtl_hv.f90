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
      !!
      !! It uses the temperature function:
      !!
      !! \[ \tau_{ij}(T) =  \frac{g_{ij}^0}{RT} + \frac{g_{ij}^T T}{RT} \]
      !! # Examples
      !!
      !! # References
      real(pr), allocatable :: b(:) 
         !! Covolume parameter
      real(pr), allocatable :: alpha(:, :) 
         !! \( \alpha \) matrix
      real(pr), allocatable :: gji0(:, :) 
         !! \( g_{ij} \) matrix
      real(pr), allocatable :: gjiT(:, :) 
         !! \( g_{ij} \) matrix
   contains
      procedure :: excess_gibbs => excess_gibbs
   end type NRTLHV

contains

   subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Calculate Excess Gibbs and its derivatives.
      use yaeos__models_ge_base, only: nrtl_hv_ge, nrtl_hv_tdep, nrtl_hv_tdep_linear
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

      call nrtl_hv_tdep_linear(T, self%gjiT/R, self%gji0/R, tau, dtaudt, dtaudt2)
      call nrtl_hv_ge(n=n, T=T,&
         b=self%b, alpha=self%alpha, &
         tau=tau, dtaudt=dtaudt, dtaudt2=dtaudt2, &
         Ge=Ge, Gen=Gen, GeT=GeT, GeT2=GeT2, GeTn=GeTn, Gen2=Gen2)
   end subroutine excess_gibbs

end module yaeos__models_ge_nrtlhv
