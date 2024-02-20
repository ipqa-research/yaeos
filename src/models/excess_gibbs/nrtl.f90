module yaoes_models_ge_nrtl
   use yaeos_constants, only: pr, R
   use yaeos_models_ge, only: GeModel

   type, extends(GeModel) :: NRTL
      real(pr), allocatable :: G(:, :)
      real(pr), allocatable :: tau(:, :)
   end type

contains
   subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Excess Gibbs and derivs procedure
      class(NRTL), intent(in) :: self !! Model
      real(pr), intent(in) ::n(:) !! Moles vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(out) :: Ge !! Excess Gibbs
      real(pr), intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), intent(out) :: Gen(size(n))
      real(pr), intent(out) :: GeTn(size(n))
      real(pr), intent(out) :: Gen2(size(n), size(n))

      real(pr) :: x(size(n)), ge_term(size(n))
      real(pr) :: xG(size(n), size(n), size(n)) !! x*G_ij matrix
      real(pr) :: xtG(size(n), size(n), size(n)) !! x*tau*G_ij matrix

      integer :: i, j, k, nc

      ! Obtain molar fractions
      x = n/sum(n)

      nc = size(n)

      Ge = 0
      ge_term = 0
      associate(G => self%G, tau=> self%tau)

         do concurrent(i=1:nc, j=1:nc, k=1:nc)
            xG(i, j, k) = x(i) * G(j, k)
            xtG(i, j, k) = x(i) * tau(j, k) * G(j, k)
         end do

         do i=1,nc
            ge_term(i) = sum(tau(:, i) * G(:, i)* x(:)) / sum(G(:, i) * x(:))
         end do

         ge = sum(x * ge_term) * R * T
         
         
         do i=1,nc
            Gen(i) = ge_term(i)
            
            do j=1,nc
                ! gen(i) = gen(i) + 
            end do
         end do

      end associate
   end subroutine
end module
