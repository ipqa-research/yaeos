module yaeos__models_ge_nrtl_hv
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel

   type, extends(GeModel) :: NRTLHV
      !! NRTL-HV model
      !! 
      !! Huron-Vidal modified NRTL model
      !! This model incldues the Cubic EoS repulsive parameters 
      !! of pure compounds.
      !!
      !! \[
      !! \sum_i n_i 
      !!  \frac
      !!  {\sum_j n_j b_j exp \left( -\alpha_{ji} \Delta U_{ji}/RT\right)\Delta U_{ji}}
      !!  {\sum_j n_j b_j exp(-\alpha_{ji} \Delta U_{ji}/RT)}
      !! \]
      real(pr), allocatable :: aij(:, :) !! \(\alpha_{ij} parameter\)
      real(pr), allocatable :: dUij(:, :) !! \(\Delta U_{ji}\)
   contains
      procedure :: excess_gibbs
   end type

contains

      subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
         !! Excess Gibbs and derivs procedure
         class(NRTLHV), intent(in) :: self !! Model
         real(pr), intent(in) ::n(:) !! Moles vector
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr), optional, intent(out) :: Ge !! Excess Gibbs
         real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
         real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
         real(pr), optional, intent(out) :: Gen(size(n))
         real(pr), optional, intent(out) :: GeTn(size(n))
         real(pr), optional, intent(out) :: Gen2(size(n), size(n))

         integer :: i, j

         real(pr) :: up, down

         real(pr) :: b(size(n))
         real(pr) :: dU(size(n), size(n))
         real(pr) :: alpha(size(n), size(n))
         real(pr) :: aux


         Ge = 0

         b = self%b
         dU = self%dUij
         alpha = self%aij

         do i=1,nc
            up = 0
            down = 0
            do j=1,nc
            end do
         end do

      end subroutine
end module