module yaeos_models_ge_group_contribution_unifac
   use yaeos_constants, only: pr, R
   use yaeos_models_ge, only: GeModel

   type, extends(GeModel) :: UNIFAC
      real(pr), allocatable :: q(:)
      real(pr), allocatable :: r(:)
      real(pr), allocatable :: z(:)
   end type

contains
   subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Excess Gibbs and derivs procedure
      class(UNIFAC), intent(in) :: self !! Model
      real(pr), intent(in) ::n(:) !! Moles vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(out) :: Ge !! Excess Gibbs
      real(pr), intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), intent(out) :: Gen(size(n))
      real(pr), intent(out) :: GeTn(size(n))
      real(pr), intent(out) :: Gen2(size(n), size(n))

      real(pr) :: x(size(n))
      real(pr) :: ln_gamma_c(size(n)), ln_gamma_r(size(n))

      real(pr) :: theta(size(n)), phi(size(n)), L(size(n))

      integer :: i, j, nc

      x = n/sum(n)

      associate(&
         q => self%q, r => self%r, z => self%z &
         )

         theta = x * q / sum(x * q)
         phi = x * r / sum(x * r)

         L = 0.5_pr * z * (r - q) - (r - 1)

         ln_gamma_c = log(phi/x) + z/2*q * ln(theta/phi) + L - phi/x * sum(x*L)
      end associate
   end subroutine
end module
