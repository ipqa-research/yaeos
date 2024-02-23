module yaeos_models_ge_group_contribution_unifac
   ! use yaeos_constants, only: pr, R
   ! use yaeos_models_ge, only: GeModel
   integer, parameter :: pr=8
   real(pr), parameter :: R=0.0831

   integer, parameter :: ng = 3

   ! Test groups:
   ! CH3, CH2, OH

   type :: Groups
      integer, allocatable :: groups_ids(:) !! Indexes (ids) of each group
      integer, allocatable :: number_of_groups(:) !! \(\nu\)
      real(pr) :: surface_area !! q
      real(pr) :: volume !! r
   end type

   type :: UNIFAC
      ! UNIFAC parameters
      integer :: ngroups !! Total number of individual groups
      real(pr) :: z = 10
      real(pr), allocatable :: Eij(:, :) !! Groups interactions

      real(pr), allocatable :: group_area(:) !! Q_k
      real(pr), allocatable :: group_volume(:) !! R_k

      type(Groups), allocatable :: molecules(:)
      type(Groups) :: groups_stew
   end type

   type :: TemperatureFunction
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
         associate (q => self%molecules%surface_area, r => self%molecules%volume, z => self%z)
            theta = x * q / sum(x * q)
            phi = x * r / sum(x * r)
            L = 0.5_pr * z * (r - q) - (r - 1)
            ln_gamma_c = log(phi/x) + z/2*q * log(theta/phi) + L - phi/x * sum(x*L)
         end associate
      end block combinatorial

      residual: block

         do i=1,size(self%molecules)

            ! Theta_i
            !TODO: This should be in the init
         end do
      end block residual
   end subroutine

   subroutine temperature_dependence(self, T, psi, dpsidt, dpsidt2)
      class(UNIFAC) :: self
      real(pr), intent(in) :: T
      real(pr), intent(out) :: psi(:, :)
      real(pr), optional, intent(out) :: dpsidt(:, :)
      real(pr), optional, intent(out) :: dpsidt2(:, :)

      integer :: i, j
      integer :: ig, jg


      do concurrent(i=1:self%ngroups, j=1:self%ngroups)
         ig = self%groups_stew%groups_ids(i)
         jg = self%groups_stew%groups_ids(j)
         print *, ig, jg, self%Eij(ig, jg)
         psi(i, j) = exp(-self%Eij(ig, jg) / T)
      end do

   end subroutine

   subroutine group_fraction(self, x, gf, dgfdx, dgfdx2)
      type(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: gf(:)
      real(pr), optional, intent(out) :: dgfdx(self%ngroups, self%ngroups)
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
               k = molecules(i)%groups_ids(j)
               gf(k) = gf(k) + x(i) * molecules(i)%number_of_groups(i)
            end do
         end do

         gf = gf/total_groups
      end associate
   end subroutine

   subroutine group_area_fraction(&
      self, x, &
      gf, dgfdx, dgfdx2, &
      theta, dthetadx, dthetadx2&
      )
      type(UNIFAC) :: self
      real(pr), intent(in) :: x(:)

      real(pr), intent(in) :: gf(:)
      real(pr), optional, intent(in) :: dgfdx(ng, ng)
      real(pr), optional, intent(in) :: dgfdx2(ng, ng, ng)

      real(pr), intent(out) :: theta(ng)
      real(pr), optional, intent(out) :: dthetadx(ng, ng)
      real(pr), optional, intent(out) :: dthetadx2(ng, ng, ng)


      real(pr) :: total_groups_area
      integer :: i, gi

      associate(group_area => self%group_area, stew => self%groups_stew)

         theta = gf * group_area

         do i=1,size(stew%groups_ids)
            gi = stew%groups_ids(i)
            total_groups_area = total_groups_area + group_area(gi) * gf(gi)
         end do

         theta = theta/total_groups_area

      end associate
   end subroutine

   subroutine group_residual()
   end subroutine

   type(UNIFAC) function setup_unifac(molecules, Eij, Qk, Rk)
      type(Groups), intent(in) :: molecules(:)
      
      real(pr), optional, intent(in) :: Eij(:, :)
      real(pr), optional, intent(in) :: Qk(:)
      real(pr), optional, intent(in) :: Rk(:)
      type(Groups) :: soup

      integer :: gi, i, j

      setup_unifac%molecules = molecules

      allocate(soup%groups_ids(0))
      allocate(soup%number_of_groups(0))

      ! Get all the groups indexes and counts into a single stew of groups.
      do i=1,size(molecules)
         do j=1,size(molecules(i)%groups_ids)
            gi = molecules(i)%groups_ids(j)

            if (all(soup%groups_ids - gi  /= 0)) then
               ! Add group if it wasn't included yet
               soup%groups_ids = [soup%groups_ids, gi]
               soup%number_of_groups = [soup%number_of_groups, 0]
            end if

            gi = findloc(soup%groups_ids - gi, 0, dim=1)

            soup%number_of_groups(gi) = soup%number_of_groups(gi) &
               + molecules(i)%number_of_groups(gi)
         end do
      end do

      setup_unifac%groups_stew = soup
      setup_unifac%ngroups = size(soup%number_of_groups)

      if (present(Eij)) setup_unifac%Eij = Eij
      if (present(Qk)) setup_unifac%group_area = Qk
      if (present(Rk)) setup_unifac%group_volume = Rk
   end function
end module

program main
   use yaeos_models_ge_group_contribution_unifac
   use stdlib_io_npy, only: load_npy
   implicit none

   type(Groups) :: molecules(2)
   type(UNIFAC) :: model
   real(pr) :: x(2) = [0.3, 0.7]
   real(pr) :: gf(ng), theta(ng)
   
   real(pr), allocatable :: Aij(:, :)
   real(pr) :: Qk(179), Rk(179)
   real(pr) :: psi(3, 3)

   integer :: i, j, gi
   call load_npy("test/data/aij.npy", Aij)
   print *, aij(14, 14)

   ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]

   ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   Qk(1) = 0.848
   Qk(2) = 0.540
   Qk(14) = 1.2

   Rk(1) = 0.9011_pr
   Rk(2) = 0.6744_pr
   Rk(14) =  1.0_pr
   model = setup_unifac(&
             molecules, &
             Eij=Aij, &
             Qk=Qk, &
             Rk=Rk &
        )

   ! call group_fraction(model, x, gf)
   ! call group_area_fraction(model, x, gf=gf, theta=theta)


   ! print *, gf
   ! print *, theta
   psi = 0
   call temperature_dependence(model, 200._pr, psi)
   print *, psi(1, :)
   print *, psi(2, :)
   print *, psi(3, :)
end program