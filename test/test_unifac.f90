module yaeos_models_ge_group_contribution_unifac
   ! use yaeos_constants, only: pr, R
   ! use yaeos_models_ge, only: GeModel
   implicit none
   integer, parameter :: pr=8
   real(pr), parameter :: R=0.0831

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
         !! Substances present in the system
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
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs
      real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), optional, intent(out) :: Gen(size(n))
      real(pr), optional, intent(out) :: GeTn(size(n))
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))

      real(pr) :: x(size(n))
      real(pr) :: ln_gamma_c(size(n)), ln_gamma_r(size(n)), ln_activity(size(n))


      integer :: i, j, nc

      x = n/sum(n)

      nc = size(self%molecules)

      ln_activity = 0

      call combinatorial_activity(self, x, ln_gamma_c)
      call residual_activity(self, x, T, ln_gamma_r)
      ln_activity = ln_gamma_c + ln_gamma_r

      print *, exp(ln_activity)
      
      if (present(Ge)) Ge = sum(x * ln_activity)
   end subroutine

   subroutine residual_activity(self, x, T, ln_gamma_r)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: ln_gamma_r(:)

      real(pr) :: ln_Gamma(size(self%groups_stew%groups_ids))
      real(pr) :: ln_Gamma_i(size(self%groups_stew%groups_ids))

      real(pr) :: xpure(size(x))
      integer :: i, k, nc

      ln_gamma = 0
      nc = size(self%molecules)

      call group_big_gamma(self, x, T, ln_Gamma)
      do i=1,nc
         xpure = 0
         xpure(i) = 1

         ln_gamma_i = 0
         call group_big_gamma(self, xpure, T, ln_Gamma_i)
         do k=1,size(self%molecules(i)%groups_ids)
            ln_gamma_r(i) = ln_gamma_r(i) + self%molecules(i)%number_of_groups(k) * (ln_gamma(k) - ln_gamma_i(k))
         end do
      end do
   end subroutine

   subroutine combinatorial_activity(self, x, ln_gamma_c)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: ln_gamma_c(:)

      real(pr) :: theta(size(x)), phi(size(x)), L(size(x))

      associate (&
         q => self%molecules%surface_area,&
         r => self%molecules%volume,&
         z => self%z &
         )

         theta = x * q / sum(x * q)
         phi = x * r / sum(x * r)
         L = 0.5_pr * z * (r - q) - (r - 1)
         ln_gamma_c = log(phi/x) + z/2*q * log(theta/phi) + L - phi/x * sum(x*L)
      end associate
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
         psi(i, j) = exp(-self%Eij(ig, jg) / T)
      end do

   end subroutine

   subroutine group_fraction(self, x, gf, dgfdx, dgfdx2)
      type(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: gf(:)
      real(pr), optional, intent(out) :: dgfdx(:, :)
      real(pr), optional, intent(out) :: dgfdx2(:, :, :)

      integer :: i, j, k, nc

      real(pr) :: total_groups, total_group_k(size(gf))
      associate(molecules => self%molecules)

         nc = size(x)

         total_group_k = self%groups_stew%number_of_groups

         total_groups = 0
         do i=1,nc
            do j=1,size(molecules(i)%number_of_groups)
               total_groups = total_groups &
                  + x(i) * molecules(i)%number_of_groups(j)
            end do
         end do

         gf = 0

         do i=1,nc
            do j=1,size(molecules(i)%number_of_groups)
               k = findloc(&
                  self%groups_stew%groups_ids,&
                  molecules(i)%groups_ids(j), dim=1 &
                  )
               
               gf(k) = gf(k) + x(i) * molecules(i)%number_of_groups(j)
               if (present(dgfdx)) then
                  dgfdx(k, i) = dgfdx(k, i) &
                              + (1 * total_group_k(k))/total_groups &
                              * molecules(i)%number_of_groups(j)
               end if
            end do
         end do

         gf = gf/total_groups
         if (present(dgfdx)) dgfdx = dgfdx/total_groups
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
      real(pr), optional, intent(in) :: dgfdx(:, :)
      real(pr), optional, intent(in) :: dgfdx2(:, :, :)

      real(pr), intent(out) :: theta(:)
      real(pr), optional, intent(out) :: dthetadx(:, :)
      real(pr), optional, intent(out) :: dthetadx2(:, :, :)


      real(pr) :: total_groups_area
      integer :: i, gi

      associate(group_area => self%group_area, stew => self%groups_stew)

         theta = gf * group_area(self%groups_stew%groups_ids)

         do i=1,size(stew%groups_ids)
            gi = stew%groups_ids(i)
            total_groups_area = total_groups_area + group_area(gi) * gf(i)
         end do

         theta = theta/total_groups_area

      end associate
   end subroutine

   subroutine group_big_gamma(self, x, T, ln_Gamma)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: ln_Gamma(:)

      real(pr) :: gf(size(self%groups_stew%groups_ids))
      real(pr) :: theta(size(self%groups_stew%groups_ids))
      real(pr) :: psi(&
         size(self%groups_stew%groups_ids),&
         size(self%groups_stew%groups_ids) &
         )

      integer :: ng
      integer :: n, m, k, gi

      real(pr) :: up, down

      ng = size(self%groups_stew%groups_ids)
      call group_fraction(self, x, gf=gf)
      call group_area_fraction(self, x, gf=gf, theta=theta)
      call temperature_dependence(self, T, psi)

      do k=1,ng

         gi = self%groups_stew%groups_ids(k)

         up = 0
         do  m=1,ng
            down = 0
            do n=1,ng
               down = down + theta(n) * psi(n, m)
            end do
            up = up + theta(m) * psi(k, m) / down
         end do

         ln_Gamma(k) = self%group_area(gi) * (&
            1 - log(sum(theta * psi(:, k))) - up &
            )

      end do
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

      associate(&
         r => setup_unifac%molecules%volume, &
         q => setup_unifac%molecules%surface_area &
         )

         ! Get all the groups indexes and counts into a single stew of groups.
         do i=1,size(molecules)

            r(i) = 0
            q(i) = 0

            do j=1,size(molecules(i)%groups_ids)
               gi = molecules(i)%groups_ids(j)

               r(i) = r(i) + molecules(i)%number_of_groups(j) * Rk(gi)
               q(i) = q(i) + molecules(i)%number_of_groups(j) * Qk(gi)

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

      end associate

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

   real(pr), allocatable :: Aij(:, :)
   real(pr), allocatable :: Qk(:), Rk(:)

   real(pr) :: psi(3, 3)
   real(pr) :: gf(3), dgdx(3, 2)

   integer :: i, j, gi

   call load_npy("data/unifac_aij.npy", Aij)
   call load_npy("data/unifac_Qk.npy", Qk)
   call load_npy("data/unifac_Rk.npy", Rk)

   ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]

   ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   model = setup_unifac(&
      molecules, &
      Eij=Aij, &
      Qk=Qk, &
      Rk=Rk &
      )

   psi = 0
   call temperature_dependence(model, 200._pr, psi)
   call excess_gibbs(model, x, 200._pr)

   dgdx = 0
   call group_fraction(model, x, gf, dgdx)

   print *, gf
   print *, dgdx(1, :)
   print *, dgdx(2, :)
   print *, dgdx(3, :)
end program
