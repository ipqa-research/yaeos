module yaeos_models_ge_group_contribution_unifac
   ! use yaeos_constants, only: pr, R
   ! use yaeos_models_ge, only: GeModel
   integer, parameter :: pr=8
   real(pr), parameter :: R=0.0831

   integer, parameter :: ng = 3

   ! Test groups:
   ! CH3, CH2, OH

   ! UNIFAC parameters
   real(pr) :: Eij(ng, ng) !! Groups interactions
   real(pr) :: group_area(ng) =   [0.848, 0.540, 1.0] !! Q_k
   real(pr) :: group_volume(ng) = [0.9011, 0.6744, 1.2] !! R_k
   type :: Groups
      integer, allocatable :: groups_indexes(:)
      integer, allocatable :: number_of_groups(:) !! \(\nu\)
      real(pr) :: surface_area = 0
      real(pr) :: volume = 0
   end type
   type :: UNIFAC
      real(pr) :: z = 10
      integer :: total_groups
      type(Groups), allocatable :: molecules(:)
   end type

contains

   subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Excess Gibbs and derivs procedure
      class(UNIFAC) :: self !! Model
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


      integer :: i, j, nc

      x = n/sum(n)

      nc = size(self%molecules)

      combinatorial: block
         real(pr) :: theta(size(n)), phi(size(n)), L(size(n))
         associate (q => self%molecules%surface_area, r => self%molecules%volume)
            theta = x * q / sum(x * q)
            phi = x * r / sum(x * r)
            L = 0.5_pr * z * (r - q) - (r - 1)
            ln_gamma_c = log(phi/x) + z/2*q * ln(theta/phi) + L - phi/x * sum(x*L)
         end associate
      end block combinatorial

      residual: block

         do i=1,size(self%molecules)
            
            ! Theta_i 
            !TODO: This should be in the init
         end do
      end block residual
   end subroutine


   subroutine group_fraction(self, x, gf, dgfdx, dgfdx2)
      type(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: gf(ng)
      real(pr), optional, intent(out) :: dgfdx(ng, ng)
      real(pr), optional, intent(out) :: dgfdx2(ng, ng, ng)

      integer :: i, j, k, nc

      real(pr) :: total_groups
      associate(molecules => self%molecules)

      nc = size(x)

      total_groups = 0
      do i=1,nc
         do j=1,size(molecules(i)%number_of_groups)
            total_groups = total_groups + x(i) * molecules(i)%number_of_groups(j)
         end do
      end do

      gf = 0
      do i=1,nc
         do j=1,size(molecules(i)%number_of_groups)
            k = molecules(i)%groups_indexes(j)
            gf(k) = gf(k) + x(i) * molecules(i)%number_of_groups(i)
         end do
      end do
      gf = gf/total_groups
      end associate
   end subroutine
end module

program main
   use yaeos_models_ge_group_contribution_unifac

   type(Groups) :: subs(2)
   type(UNIFAC) :: model


   subs(1)%groups_indexes = [1]
   subs(1)%number_of_groups = [2]




end program
