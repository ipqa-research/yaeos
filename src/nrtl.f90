module yaeos__models_ge_NRTL
   !! # NRTL activity-coefficient model
   !!
   !! # Description
   !! Non-Random-Two-Liquid model for excess Gibbs energy. The constructor is
   !! [[NRTL]] and the model extends [[GeModel]]. It uses local compositions
   !! and temperature-dependent pair interactions to represent liquid
   !! non-ideality.
   !!
   !! # Examples
   !! ```fortran
   !! model = NRTL(A, B, alpha)
   !! call model%ln_activity_coefficient(n, T, ln_gamma)
   !! ```
   !!
   !! # References
   !! Renon, H., and J. M. Prausnitz (1968). Local compositions in
   !! thermodynamic excess functions for liquid mixtures. *AIChE Journal*, 14,
   !! 135-144.
   use yaeos__models_ge, only: GeModel
   use yaeos__constants, only: pr, R
   implicit none

   type, extends(GeModel) ::  NRTL
      !! # NRTL
      !! Non-Random-Two-Liquid excess-Gibbs model.
      !!
      !! # Description
      !! \[
      !!    G^E = nRT \cdot \sum_i x_i \frac{\sum_j x_j \tau_{ji} G_{ji}}{\sum_j x_j G_{ji}}
      !! \]
      !!
      !! with:
      !!
      !! \[\tau_{ij} = A_{ij} + \frac{B_{ij}}{T}\]
      !!
      !! \[G_{ij} = \exp(-C_{ij}\tau_{ij})\]
      !!
      !! # Examples
      !! ```fortran
      !! model = NRTL(A, B, C)
      !! ```
      real(pr), allocatable :: a(:, :) !! \(A_{ij}\) matrix
      real(pr), allocatable :: b(:, :) !! \(B_{ij}\) matrix [K]
      real(pr), allocatable :: c(:, :) !! \(C_{ij}\) non-randomness matrix
   contains
      procedure :: excess_gibbs => excess_gibbs
   end type NRTL

   interface NRTL
      module procedure :: init
   end interface NRTL
contains
   type(NRTL) function init(a, b, c)
      !! # init
      !! Construct an NRTL model from its interaction matrices.
      !!
      !! # Description
      !! Stores the \(A_{ij}\) (dimensionless), \(B_{ij}\) [K], and
      !! non-randomness \(C_{ij}\) matrices used by the NRTL temperature and
      !! local-composition terms.
      !!
      !! # Examples
      !! ```fortran
      !! model = NRTL(A, B, C)
      !! ```
      real(pr), intent(in) :: a(:, :) !! \(A_{ij}\) matrix
      real(pr), intent(in) :: b(:, :) !! \(B_{ij}\) matrix
      real(pr), intent(in) :: c(:, :) !! \(C_{ij}\) matrix

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
      real(pr), optional, intent(out) :: GeTn(size(n)) !! Mixed temperature-mole-number derivative.
      real(pr), optional, intent(out) :: Gen2(size(n), size(n)) !! Hessian with respect to mole amounts.
      real(pr) :: tau(size(n), size(n))
      real(pr) :: dtaudt(size(n), size(n))
      real(pr) :: dtaudt2(size(n), size(n))

      real(pr) :: b(size(n))

      b = 1

      call nrtl_hv_tdep_linear(T, self%a, self%b, tau, dtaudt, dtaudt2)
      call nrtl_hv_ge(n=n, T=T,&
         b=b, alpha=self%c, &
         tau=tau, dtaudt=dtaudt, dtaudt2=dtaudt2, &
         Ge=Ge, Gen=Gen, GeT=GeT, GeT2=GeT2, GeTn=GeTn, Gen2=Gen2)
   end subroutine excess_gibbs

end module yaeos__models_ge_NRTL
