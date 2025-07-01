module yaeos__models_ge_NRTL
   use yaeos__models_ge, only: GeModel
   use yaeos__tapenade_interfaces
   use yaeos__constants, only: pr, R
   implicit none

   type, extends(GeModel) ::  NRTL
      !! Non-Random-Two-Liquid model
      !!
      !! \[
      !!    G^E = nRT \cdot \sum_i x_i \frac{\sum_j x_j \tau_{ji} G_{ji}}{\sum_j x_j G_{ji}}
      !! \]
      !!
      !! with:
      !!
      !! \[\tau_{ij} = A_{ij} + \frac{B_{ij}}{T}\]
      !!
      !! \[G_{ij} = exp(-\frac{C_{ij}}{\tau_{ij}})\]
      real(pr), allocatable :: a(:, :) !! A_{ij} matrix
      real(pr), allocatable :: b(:, :) !! B_{ij} matrix
      real(pr), allocatable :: c(:, :) !! C_{ij} matrix
   contains
      procedure :: excess_gibbs => excess_gibbs
   end type NRTL

   interface NRTL
      module procedure :: init
   end interface NRTL
contains
   type(NRTL) function init(a, b, c)
      real(pr), intent(in) :: a(:, :)
      real(pr), intent(in) :: b(:, :)
      real(pr), intent(in) :: c(:, :)

      init%a = a
      init%b = b
      init%c = c
   end function init

   subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Calculate Excess Gibbs and its derivatives.
      use yaeos__models_ge_base, only: nrtl_hv_ge, nrtl_hv_tdep_linear
      class(NRTL), intent(in) :: self !! Model
      real(pr), intent(in) ::n(:) !! Moles vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs free energy
      real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), optional, intent(out) :: Gen(size(n)) !! \(\frac{dG^E}{dn}\)
      real(pr), optional, intent(out) :: GeTn(size(n))
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))
      real(pr) :: tau(size(n), size(n)), dtaudt(size(n), size(n)), dtaudt2(size(n), size(n))

      real(pr) :: b(size(n))

      b = 1

      call nrtl_hv_tdep_linear(T, self%a, self%b, tau, dtaudt, dtaudt2)
      call nrtl_hv_ge(n=n, T=T,&
         b=b, alpha=self%c, &
         tau=tau, dtaudt=dtaudt, dtaudt2=dtaudt2, &
         Ge=Ge, Gen=Gen, GeT=GeT, GeT2=GeT2, GeTn=GeTn, Gen2=Gen2)
   end subroutine excess_gibbs

end module yaeos__models_ge_NRTL